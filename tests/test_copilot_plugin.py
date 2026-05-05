from tas_copilot_contracts import ToolCall, ToolResultStatus

from rescalculator.copilot_plugin import get_tool_specs, run_tool


def _arguments(backend="numpy"):
    return {
        "q_point": {"h": 1.0, "k": 0.0, "l": 0.0},
        "energy_transfer": {"value": {"value": 5.0, "unit": "meV"}},
        "instrument_config": {
            "fixed_energy_mode": "fixed_final_energy",
            "fixed_energy": {"value": 14.7, "unit": "meV"},
            "monochromator": "PG(002)",
            "analyzer": "PG(002)",
            "collimation": "40-40-40-80",
        },
        "lattice": {
            "a": 5.0,
            "b": 5.0,
            "c": 5.0,
            "alpha": 90.0,
            "beta": 90.0,
            "gamma": 90.0,
            "orient1": [[1, 0, 0]],
            "orient2": [[0, 1, 0]],
        },
        "backend": backend,
    }


def _call(namespace, name, arguments):
    return ToolCall(
        tool_namespace=namespace,
        tool_name=name,
        arguments=arguments,
        caller_workflow="test",
    )


def test_get_tool_specs_exposes_resolution_and_angle_tools():
    tools = {(spec.namespace, spec.name) for spec in get_tool_specs()}

    assert ("tas.resolution", "calculate_resolution") in tools
    assert ("tas.device_angles", "calculate_angles") in tools


def test_auto_backend_prefers_numba_when_installed():
    import importlib.util

    from rescalculator.backends import get_backend

    backend = get_backend("auto")
    assert backend.name in {"numba", "numpy"}
    if importlib.util.find_spec("numba") is not None:
        assert backend.name == "numba"


def test_calculate_resolution_returns_matrix_and_backend():
    result = run_tool(_call("tas.resolution", "calculate_resolution", _arguments()))

    assert result.status == ToolResultStatus.success
    assert result.output["backend_used"] == "numpy"
    assert result.output["prefactor"] > 0
    matrix = result.output["resolution"]["resolution_matrix"]
    assert len(matrix) == 4
    assert len(matrix[0]) == 4


def test_calculate_angles_returns_motor_positions():
    result = run_tool(_call("tas.device_angles", "calculate_angles", _arguments()))

    assert result.status == ToolResultStatus.success
    angles = result.output["device_angles"]["angles"]
    assert [item["name"] for item in angles] == ["M1", "M2", "S1", "S2", "A1", "A2"]
    assert all(item["position"]["unit"] == "deg" for item in angles)


def test_unknown_backend_reports_unavailable_or_error():
    args = _arguments(backend="missing_backend")
    result = run_tool(_call("tas.resolution", "calculate_resolution", args))

    assert result.status == ToolResultStatus.error
