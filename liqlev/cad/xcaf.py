from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path

import cadquery as cq
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.TCollection import TCollection_AsciiString, TCollection_ExtendedString
from OCP.TDataStd import TDataStd_Name
from OCP.TDF import TDF_Label, TDF_LabelSequence, TDF_Tool
from OCP.TDocStd import TDocStd_Document
from OCP.TopLoc import TopLoc_Location
from OCP.XCAFDoc import XCAFDoc_DocumentTool, XCAFDoc_ShapeTool


class StepProductError(RuntimeError):
    """Raised when an exact, usable STEP product cannot be selected."""


@dataclass(frozen=True)
class _ProductMatch:
    matched_label_entry: str
    instance_label_entry: str
    product_name: str
    shape: cq.Shape


def sha256_file(path: str | Path) -> str:
    """Return the file's SHA-256 digest as uppercase hexadecimal."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as file_obj:
        for block in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _label_entry(label: TDF_Label) -> str:
    entry = TCollection_AsciiString()
    TDF_Tool.Entry_s(label, entry)
    return entry.ToCString()


def _label_name(label: TDF_Label) -> str | None:
    name_attribute = TDataStd_Name()
    if not label.FindAttribute(TDataStd_Name.GetID_s(), name_attribute):
        return None
    return name_attribute.Get().ToExtString()


def _resolved_label(label: TDF_Label) -> TDF_Label:
    referred_label = TDF_Label()
    if XCAFDoc_ShapeTool.GetReferredShape_s(label, referred_label):
        return referred_label
    return label


def _traverse_named_products(
    shape_tool: XCAFDoc_ShapeTool,
    product_name: str,
) -> tuple[_ProductMatch, ...]:
    matches: list[_ProductMatch] = []

    def visit(label: TDF_Label, parent_location: TopLoc_Location) -> None:
        local_location = XCAFDoc_ShapeTool.GetLocation_s(label)
        complete_location = parent_location.Multiplied(local_location)
        resolved_label = _resolved_label(label)
        resolved_name = _label_name(resolved_label)

        if resolved_name == product_name:
            resolved_shape = XCAFDoc_ShapeTool.GetShape_s(resolved_label)
            if resolved_shape.IsNull():
                raise StepProductError(
                    f"resolved product {product_name!r} has a null shape"
                )
            placed_shape = cq.Shape.cast(resolved_shape.Located(complete_location))
            if placed_shape.isNull():
                raise StepProductError(
                    f"resolved product {product_name!r} has an unusable shape"
                )
            matches.append(
                _ProductMatch(
                    matched_label_entry=_label_entry(resolved_label),
                    instance_label_entry=_label_entry(label),
                    product_name=resolved_name,
                    shape=placed_shape,
                )
            )

        components = TDF_LabelSequence()
        if XCAFDoc_ShapeTool.GetComponents_s(
            resolved_label,
            components,
            False,
        ):
            for index in range(1, components.Length() + 1):
                visit(components.Value(index), complete_location)

    free_shapes = TDF_LabelSequence()
    shape_tool.GetFreeShapes(free_shapes)
    for index in range(1, free_shapes.Length() + 1):
        visit(free_shapes.Value(index), TopLoc_Location())
    return tuple(matches)


def _select_unique_match(
    matches: tuple[_ProductMatch, ...],
    *,
    product_name: str,
    source: Path,
) -> _ProductMatch:
    if not matches:
        raise StepProductError(
            f"zero exact product matches for {product_name!r} in {source}"
        )
    if len(matches) != 1:
        entries = ", ".join(match.instance_label_entry for match in matches)
        raise StepProductError(
            f"multiple exact product matches ({len(matches)}) for "
            f"{product_name!r} in {source}: {entries}"
        )
    return matches[0]


def load_named_product(
    path: str | Path,
    *,
    product_name: str,
    expected_sha256: str,
) -> cq.Shape:
    """Load exactly one named XCAF product with its assembly placement."""

    try:
        source = Path(path)
        actual_sha256 = sha256_file(source)
    except (OSError, TypeError, ValueError) as exc:
        raise StepProductError(f"STEP source is missing or unreadable: {path}") from exc

    if actual_sha256 != expected_sha256.upper():
        raise StepProductError(
            "STEP source SHA-256 mismatch: "
            f"expected {expected_sha256.upper()}, got {actual_sha256}"
        )

    try:
        document = TDocStd_Document(TCollection_ExtendedString("MDTV-XCAF"))
        reader = STEPCAFControl_Reader()
        reader.SetNameMode(True)
        read_status = reader.ReadFile(str(source))
        if read_status != IFSelect_RetDone:
            raise StepProductError(
                f"STEP read failed for {source}: status {read_status}"
            )
        if not reader.Transfer(document):
            raise StepProductError(f"STEP transfer failed for {source}")
        shape_tool = XCAFDoc_DocumentTool.ShapeTool_s(document.Main())
        matches = _traverse_named_products(shape_tool, product_name)
    except StepProductError:
        raise
    except Exception as exc:
        raise StepProductError(
            f"STEP read, transfer, or XCAF traversal failed for {source}"
        ) from exc

    match = _select_unique_match(
        matches,
        product_name=product_name,
        source=source,
    )
    return match.shape
