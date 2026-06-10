"""TAS co-pilot plugin entry points for rescalculator."""

from __future__ import annotations

from typing import Any

import numpy as np


def get_tool_specs():
    from tas_copilot_contracts import Availability, ImplementationType, PermissionLevel, SideEffectLevel, ToolSpec

    common = {
        "version": "2.0.1",
        "availability": Availability.available,
        "implementation_type": ImplementationType.python_entrypoint,
        "required_dependencies": ["rescalculator", "icp-lattice-calculator"],
        "provenance_policy": "Record Q/E/lattice/instrument configuration, backend, rescalculator version, and output hashes.",
        "enabled_by_default": True,
        "owner_repo": "git@github.com:scattering/rescalculator.git",
        "documentation_url": "https://github.com/scattering/rescalculator",
        "maintainer": "William Ratcliff",
    }
    return [
        ToolSpec(
            name="calculate_resolution",
            namespace="tas.resolution",
            description="Calculate a TAS resolution matrix using the external rescalculator package.",
            input_schema={
                "type": "object",
                "required": ["q_point", "energy_transfer", "instrument_config", "lattice"],
                "properties": {
                    "q_point": {"type": "object"},
                    "energy_transfer": {"type": "object"},
                    "instrument_config": {"type": "object"},
                    "lattice": {"type": "object"},
                    "backend": {"type": "string"},
                    "experiment": {"type": "object"},
                },
            },
            output_schema={"type": "object"},
            side_effect_level=SideEffectLevel.none,
            preconditions=["Q/E, lattice, orientation, and spectrometer configuration are supplied."],
            postconditions=["No files or instrument state are modified."],
            permissions_required=[PermissionLevel.read_only],
            **common,
        ),
        ToolSpec(
            name="calculate_angles",
            namespace="tas.device_angles",
            description="Calculate TAS motor angles for a requested Q/E using icp-lattice-calculator.",
            input_schema={
                "type": "object",
                "required": ["q_point", "energy_transfer", "instrument_config", "lattice"],
                "properties": {
                    "q_point": {"type": "object"},
                    "energy_transfer": {"type": "object"},
                    "instrument_config": {"type": "object"},
                    "lattice": {"type": "object"},
                    "experiment": {"type": "object"},
                },
            },
            output_schema={"type": "object"},
            side_effect_level=SideEffectLevel.none,
            preconditions=["Q/E, lattice, orientation, and spectrometer configuration are supplied."],
            postconditions=["No files or instrument state are modified."],
            permissions_required=[PermissionLevel.read_only],
            **common,
        ),
        ToolSpec(
            name="calculate_ub",
            namespace="tas.ubmatrix",
            description="Build an orientation-aware UB-like matrix and reciprocal-lattice summary using icp-lattice-calculator.",
            input_schema={
                "type": "object",
                "required": ["lattice"],
                "properties": {
                    "lattice": {"type": "object"},
                    "ub_id": {"type": "string"},
                },
            },
            output_schema={"type": "object"},
            side_effect_level=SideEffectLevel.none,
            preconditions=["Lattice parameters and two orientation vectors are supplied."],
            postconditions=["UBMatrix-compatible output and reciprocal-lattice summary are returned."],
            permissions_required=[PermissionLevel.read_only],
            **common,
        ),
        ToolSpec(
            name="calculate_reciprocal_angle",
            namespace="tas.ubmatrix",
            description="Angle between two reciprocal-space (HKL) vectors from the icp-lattice-calculator metric; not constrained to a scattering plane.",
            input_schema={
                "type": "object",
                "required": ["lattice", "v1", "v2"],
                "properties": {
                    "lattice": {"type": "object"},
                    "v1": {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3},
                    "v2": {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3},
                },
            },
            output_schema={"type": "object"},
            side_effect_level=SideEffectLevel.none,
            preconditions=["Lattice parameters and two HKL vectors are supplied."],
            postconditions=["No files or instrument state are modified."],
            permissions_required=[PermissionLevel.read_only],
            **common,
        ),
        ToolSpec(
            name="calculate_inplane_orthogonal",
            namespace="tas.ubmatrix",
            description="In-plane perpendicular HKL direction for a vector in the (orient1, orient2) scattering plane, under the reciprocal-lattice metric (icp-lattice-calculator StandardSystem).",
            input_schema={
                "type": "object",
                "required": ["lattice", "v"],
                "properties": {
                    "lattice": {"type": "object"},
                    "v": {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3},
                },
            },
            output_schema={"type": "object"},
            side_effect_level=SideEffectLevel.none,
            preconditions=["Lattice parameters with orient1/orient2 plane vectors and one in-plane HKL vector are supplied."],
            postconditions=["No files or instrument state are modified."],
            permissions_required=[PermissionLevel.read_only],
            **common,
        ),
    ]


def run_tool(call):
    from tas_copilot_contracts import ProvenanceRecord, ToolError, ToolResult, ToolResultStatus

    dispatch = {
        ("tas.resolution", "calculate_resolution"): _calculate_resolution,
        ("tas.device_angles", "calculate_angles"): _calculate_angles,
        ("tas.ubmatrix", "calculate_ub"): _calculate_ub,
        ("tas.ubmatrix", "calculate_reciprocal_angle"): _calculate_reciprocal_angle,
        ("tas.ubmatrix", "calculate_inplane_orthogonal"): _calculate_inplane_orthogonal,
    }
    handler = dispatch.get((call.tool_namespace, call.tool_name))
    if handler is None:
        return _result(
            call,
            ToolResult,
            ToolResultStatus.validation_failed,
            errors=[ToolError(code="unknown_tool", message=f"Unknown rescalculator tool {call.tool_namespace}.{call.tool_name}")],
        )
    try:
        return handler(call, ToolResult, ToolResultStatus, ToolError, ProvenanceRecord)
    except ImportError as exc:
        return _result(call, ToolResult, ToolResultStatus.unavailable, errors=[ToolError(code="optional_dependency_missing", message=str(exc))])
    except Exception as exc:
        return _result(call, ToolResult, ToolResultStatus.error, errors=[ToolError(code=exc.__class__.__name__, message=str(exc))])


def _calculate_resolution(call, ToolResult, ToolResultStatus, ToolError, ProvenanceRecord):
    from rescalculator import TASResolution
    from tas_copilot_contracts import ResolutionResult, UnitValue

    lattice = _build_lattice(call.arguments)
    exp = _build_exp(call.arguments)
    h, k, l = _q_components(call.arguments)
    energy = _energy_transfer(call.arguments)
    backend_name = str(call.arguments.get("backend", "auto"))
    calculator = TASResolution(lattice, backend=backend_name)
    actual_backend = getattr(calculator.backend, "name", calculator.backend.__class__.__name__)
    h_arr = np.asarray([h], dtype=float)
    k_arr = np.asarray([k], dtype=float)
    l_arr = np.asarray([l], dtype=float)
    w_arr = np.asarray([energy], dtype=float)
    lattice.npts = 1
    prefactor, matrices = calculator.ResMatS(h_arr, k_arr, l_arr, w_arr, [exp])
    matrix = np.asarray(matrices[:, :, 0], dtype=float)
    diagonal = np.diag(matrix)
    fwhm = {
        f"axis_{idx}": UnitValue(value=float(2.354820045 / np.sqrt(value)), unit="rlu_or_meV")
        for idx, value in enumerate(diagonal)
        if value > 0
    }
    result = ResolutionResult(
        summary=f"Resolution calculated with rescalculator backend {actual_backend}.",
        resolution_matrix=matrix.tolist(),
        fwhm=fwhm,
        provenance=[ProvenanceRecord(source_type="tool", source_id="rescalculator", method="TASResolution.ResMatS")],
    )
    return _result(
        call,
        ToolResult,
        ToolResultStatus.success,
        output={
            "resolution": result.model_dump(mode="json"),
            "backend_requested": backend_name,
            "backend_used": actual_backend,
            "prefactor": float(np.asarray(prefactor).reshape(-1)[0]),
            "experiment": _json_ready_exp(exp),
        },
        provenance=[ProvenanceRecord(source_type="tool", source_id="rescalculator", method="TASResolution.ResMatS")],
    )


def _calculate_angles(call, ToolResult, ToolResultStatus, ToolError, ProvenanceRecord):
    from lattice_calculator import Orientation, SpecGoTo
    from tas_copilot_contracts import DeviceAngleResult, MotorPosition, UnitValue

    lattice = _build_lattice(call.arguments)
    orientation = Orientation(np.asarray(_lattice_args(call.arguments)["orient1"], dtype=float), np.asarray(_lattice_args(call.arguments)["orient2"], dtype=float))
    exp = _build_exp(call.arguments)
    h, k, l = _q_components(call.arguments)
    energy = _energy_transfer(call.arguments)
    angles = SpecGoTo(
        np.asarray([h], dtype=float),
        np.asarray([k], dtype=float),
        np.asarray([l], dtype=float),
        np.asarray([energy], dtype=float),
        [exp],
        lattice,
        orientation,
    )
    names = ["M1", "M2", "S1", "S2", "A1", "A2"]
    positions = [
        MotorPosition(name=name, position=UnitValue(value=float(np.asarray(value).reshape(-1)[0]), unit="deg"))
        for name, value in zip(names, angles)
    ]
    warnings = ["One or more calculated angles are NaN; requested Q/E may be inaccessible."] if any(np.isnan(pos.position.value) for pos in positions) else []
    result = DeviceAngleResult(
        angles=positions,
        warnings=warnings,
        provenance=[ProvenanceRecord(source_type="tool", source_id="icp-lattice-calculator", method="SpecGoTo")],
    )
    return _result(
        call,
        ToolResult,
        ToolResultStatus.success,
        output={"device_angles": result.model_dump(mode="json"), "experiment": _json_ready_exp(exp)},
        provenance=[ProvenanceRecord(source_type="tool", source_id="icp-lattice-calculator", method="SpecGoTo")],
    )


def _calculate_ub(call, ToolResult, ToolResultStatus, ToolError, ProvenanceRecord):
    from tas_copilot_contracts import UBMatrix, UnitValue

    lattice = _build_lattice(call.arguments)
    matrix = _ub_like_matrix(lattice)
    lattice_args = _lattice_args(call.arguments)
    reciprocal = {
        "astar": float(np.asarray(lattice.astar).reshape(-1)[0]),
        "bstar": float(np.asarray(lattice.bstar).reshape(-1)[0]),
        "cstar": float(np.asarray(lattice.cstar).reshape(-1)[0]),
        "alphastar": float(np.asarray(lattice.alphastar).reshape(-1)[0]),
        "betastar": float(np.asarray(lattice.betastar).reshape(-1)[0]),
        "gammastar": float(np.asarray(lattice.gammastar).reshape(-1)[0]),
    }
    result = UBMatrix(
        ub_id=str(call.arguments.get("ub_id") or "ub-from-lattice-orientation"),
        matrix=matrix.tolist(),
        lattice_parameters={
            "a": UnitValue(value=float(lattice_args["a"]), unit="angstrom"),
            "b": UnitValue(value=float(lattice_args["b"]), unit="angstrom"),
            "c": UnitValue(value=float(lattice_args["c"]), unit="angstrom"),
            "alpha": UnitValue(value=float(lattice_args["alpha"]), unit="deg"),
            "beta": UnitValue(value=float(lattice_args["beta"]), unit="deg"),
            "gamma": UnitValue(value=float(lattice_args["gamma"]), unit="deg"),
        },
        provenance=[ProvenanceRecord(source_type="tool", source_id="icp-lattice-calculator", method="Lattice reciprocal basis/orientation")],
    )
    return _result(
        call,
        ToolResult,
        ToolResultStatus.success,
        output={
            "ub_matrix": result.model_dump(mode="json"),
            "reciprocal_lattice": reciprocal,
            "orientation_basis": {
                "x": np.asarray(lattice.x, dtype=float).reshape(3).tolist(),
                "y": np.asarray(lattice.y, dtype=float).reshape(3).tolist(),
                "z": np.asarray(lattice.z, dtype=float).reshape(3).tolist(),
            },
        },
        provenance=[ProvenanceRecord(source_type="tool", source_id="icp-lattice-calculator", method="Lattice reciprocal basis/orientation")],
    )


def _calculate_reciprocal_angle(call, ToolResult, ToolResultStatus, ToolError, ProvenanceRecord):
    from lattice_calculator import angle, modvec

    arguments = call.arguments
    v1 = [float(value) for value in arguments["v1"]]
    v2 = [float(value) for value in arguments["v2"]]
    lattice = _build_lattice(arguments)
    lattice.star()
    angle_rad = _scalarize(angle(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], "latticestar", lattice))
    q1 = _scalarize(modvec(v1[0], v1[1], v1[2], "latticestar", lattice))
    q2 = _scalarize(modvec(v2[0], v2[1], v2[2], "latticestar", lattice))
    output = {
        "v1": v1,
        "v2": v2,
        "angle_deg": float(np.degrees(angle_rad)),
        "v1_modulus_inv_angstrom": q1,
        "v2_modulus_inv_angstrom": q2,
        "v1_d_spacing_angstrom": (2.0 * np.pi / q1) if q1 > 0 else None,
        "v2_d_spacing_angstrom": (2.0 * np.pi / q2) if q2 > 0 else None,
        "convention": "Reciprocal-lattice metric (latticestar) from icp-lattice-calculator; vectors are general HKL and need not lie in any scattering plane.",
    }
    return _result(
        call,
        ToolResult,
        ToolResultStatus.success,
        output=output,
        provenance=[ProvenanceRecord(source_type="tool", source_id="icp-lattice-calculator", method="lattice_calculator.angle latticestar metric")],
    )


def _calculate_inplane_orthogonal(call, ToolResult, ToolResultStatus, ToolError, ProvenanceRecord):
    from lattice_calculator import S2R, StandardSystem, angle, modvec, scalar

    arguments = call.arguments
    v = [float(value) for value in arguments["v"]]
    lattice = _build_lattice(arguments)
    lattice.star()
    lattice_args = _lattice_args(arguments)
    orient1 = np.asarray(lattice_args["orient1"], dtype=float).reshape(1, 3)
    orient2 = np.asarray(lattice_args["orient2"], dtype=float).reshape(1, 3)
    x_v, y_v, z_v = StandardSystem(orient1, orient2, lattice)
    qx = _scalarize(scalar(v[0], v[1], v[2], x_v[0, :], x_v[1, :], x_v[2, :], "latticestar", lattice))
    qy = _scalarize(scalar(v[0], v[1], v[2], y_v[0, :], y_v[1, :], y_v[2, :], "latticestar", lattice))
    qz = _scalarize(scalar(v[0], v[1], v[2], z_v[0, :], z_v[1, :], z_v[2, :], "latticestar", lattice))
    modulus = _scalarize(modvec(v[0], v[1], v[2], "latticestar", lattice))
    if modulus <= 0:
        return _result(call, ToolResult, ToolResultStatus.validation_failed, errors=[ToolError(code="zero_vector", message="v must be a nonzero HKL vector.")])
    out_of_plane_deg = float(np.degrees(np.arcsin(min(1.0, abs(qz) / modulus))))
    # Rotate the in-plane component by +90 degrees in the orthonormal plane
    # basis, then convert back to HKL through S2R. This is the metric-correct
    # in-plane perpendicular; for hexagonal cells it is NOT the naive (010).
    h_arr, k_arr, l_arr, _ = S2R(np.array([-qy]), np.array([qx]), np.array([0.0]), x_v, y_v, z_v)
    orthogonal = [_scalarize(h_arr), _scalarize(k_arr), _scalarize(l_arr)]
    rationalized = _rationalize_hkl(orthogonal)
    verification_deg = float(
        np.degrees(_scalarize(angle(v[0], v[1], v[2], orthogonal[0], orthogonal[1], orthogonal[2], "latticestar", lattice)))
    )
    warnings = []
    if out_of_plane_deg > 0.01:
        warnings.append(
            f"Input vector is {out_of_plane_deg:.3f} deg out of the (orient1, orient2) plane; the orthogonal direction was computed from its in-plane projection."
        )
    output = {
        "v": v,
        "orient1": orient1.reshape(3).tolist(),
        "orient2": orient2.reshape(3).tolist(),
        "v_out_of_plane_deg": out_of_plane_deg,
        "orthogonal_hkl": orthogonal,
        "orthogonal_hkl_rational": rationalized,
        "orthogonality_check_deg": verification_deg,
        "v_modulus_inv_angstrom": modulus,
        "orthogonal_modulus_inv_angstrom": _scalarize(modvec(orthogonal[0], orthogonal[1], orthogonal[2], "latticestar", lattice)),
        "convention": "Orthogonal within the scattering plane under the reciprocal-lattice metric (icp-lattice-calculator StandardSystem); Euclidean HKL intuition fails for non-orthogonal cells.",
        "warnings": warnings,
    }
    return _result(
        call,
        ToolResult,
        ToolResultStatus.success,
        output=output,
        provenance=[ProvenanceRecord(source_type="tool", source_id="icp-lattice-calculator", method="StandardSystem in-plane 90-degree rotation")],
    )


def _scalarize(value: Any) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def _rationalize_hkl(hkl: list[float], *, max_multiplier: int = 24, tolerance: float = 5.0e-4) -> list[int] | None:
    """Smallest integer-HKL direction parallel to the input, if one exists."""
    largest = max(abs(value) for value in hkl)
    if largest <= 0:
        return None
    normalized = [value / largest for value in hkl]
    for multiplier in range(1, max_multiplier + 1):
        scaled = [value * multiplier for value in normalized]
        rounded = [round(value) for value in scaled]
        if all(abs(s - r) <= tolerance * multiplier for s, r in zip(scaled, rounded)) and any(rounded):
            divisor = np.gcd.reduce([abs(value) for value in rounded if value]) or 1
            return [int(value / divisor) for value in rounded]
    return None


def _ub_like_matrix(lattice) -> np.ndarray:
    # Columns are the orthonormal scattering-frame basis vectors in reciprocal-lattice coordinates.
    return np.column_stack(
        [
            np.asarray(lattice.x, dtype=float).reshape(3),
            np.asarray(lattice.y, dtype=float).reshape(3),
            np.asarray(lattice.z, dtype=float).reshape(3),
        ]
    )


def _q_components(arguments: dict[str, Any]) -> tuple[float, float, float]:
    q_point = arguments.get("q_point") or {}
    return float(q_point["h"]), float(q_point["k"]), float(q_point["l"])


def _energy_transfer(arguments: dict[str, Any]) -> float:
    energy = arguments.get("energy_transfer") or {}
    value = energy.get("value", energy)
    if isinstance(value, dict):
        return float(value["value"])
    return float(value)


def _lattice_args(arguments: dict[str, Any]) -> dict[str, Any]:
    lattice = dict(arguments.get("lattice") or {})
    lattice.setdefault("a", 2 * np.pi)
    lattice.setdefault("b", 2 * np.pi)
    lattice.setdefault("c", 2 * np.pi)
    lattice.setdefault("alpha", 90.0)
    lattice.setdefault("beta", 90.0)
    lattice.setdefault("gamma", 90.0)
    lattice.setdefault("orient1", [[1, 0, 0]])
    lattice.setdefault("orient2", [[0, 1, 0]])
    return lattice


def _build_lattice(arguments: dict[str, Any]):
    from lattice_calculator import Lattice

    lattice = _lattice_args(arguments)
    return Lattice(
        a=float(lattice["a"]),
        b=float(lattice["b"]),
        c=float(lattice["c"]),
        alpha=float(lattice["alpha"]),
        beta=float(lattice["beta"]),
        gamma=float(lattice["gamma"]),
        orient1=np.asarray(lattice["orient1"], dtype=float),
        orient2=np.asarray(lattice["orient2"], dtype=float),
    )


def _build_exp(arguments: dict[str, Any]) -> dict[str, Any]:
    instrument = arguments.get("instrument_config") or {}
    overrides = arguments.get("experiment") or {}
    fixed_energy = instrument.get("fixed_energy") or {}
    fixed_energy_value = fixed_energy.get("value") if isinstance(fixed_energy, dict) else fixed_energy
    mode = str(instrument.get("fixed_energy_mode") or "fixed_final_energy").lower()
    exp = {
        "efixed": float(overrides.get("efixed", fixed_energy_value or 14.7)),
        "infin": int(overrides.get("infin", 1 if "incident" in mode else -1)),
        "dir1": int(overrides.get("dir1", 1)),
        "dir2": int(overrides.get("dir2", 1)),
        "hcol": np.asarray(overrides.get("hcol", _parse_collimation(instrument.get("collimation"))), dtype=float),
        "vcol": np.asarray(overrides.get("vcol", [120.0, 120.0, 120.0, 120.0]), dtype=float),
        "mono": {"tau": _crystal(instrument.get("monochromator"), "pg(002)"), "mosaic": float(overrides.get("mono_mosaic", 30.0))},
        "ana": {"tau": _crystal(instrument.get("analyzer"), "pg(002)"), "mosaic": float(overrides.get("ana_mosaic", 30.0))},
        "sample": {"mosaic": float(overrides.get("sample_mosaic", 30.0)), "vmosaic": float(overrides.get("sample_vmosaic", 30.0))},
        "arms": list(overrides.get("arms", [200.0, 200.0, 150.0, 150.0, 100.0])),
        "horifoc": int(overrides.get("horifoc", -1)),
        "moncor": int(overrides.get("moncor", 1)),
        "method": int(overrides.get("method", 0)),
    }
    return exp


def _parse_collimation(value: Any) -> list[float]:
    if value is None:
        return [40.0, 40.0, 40.0, 80.0]
    if isinstance(value, str):
        return [float(part.strip()) for part in value.replace("-", " ").split()[:4]]
    return [float(part) for part in value]


def _crystal(value: Any, default: str) -> str:
    if not value:
        return default
    text = str(value).strip().lower().replace(" ", "")
    aliases = {
        "pg002": "pg(002)",
        "pg(002)": "pg(002)",
        "pg004": "pg(004)",
        "pg(004)": "pg(004)",
        "ge111": "ge(111)",
        "ge(111)": "ge(111)",
        "ge220": "ge(220)",
        "ge(220)": "ge(220)",
        "ge311": "ge(311)",
        "ge(311)": "ge(311)",
        "be002": "be(002)",
        "be(002)": "be(002)",
    }
    return aliases.get(text, text)


def _json_ready_exp(exp: dict[str, Any]) -> dict[str, Any]:
    return {key: (value.tolist() if hasattr(value, "tolist") else value) for key, value in exp.items()}


def _result(call, ToolResult, status, output=None, errors=None, provenance=None):
    return ToolResult(
        status=status,
        output=output,
        errors=errors or [],
        tool_name=call.tool_name,
        tool_namespace=call.tool_namespace,
        tool_version="2.0.1",
        provenance=provenance or [],
    )
