import numpy as np
import pytest
from pathlib import Path
import pickle
import os
from neat.models.gc_bias_model import GCBiasModel, get_uniform_gc_model

def test_gc_bias_model_init():
    weights = [1.0] * 101
    model = GCBiasModel(weights, window_size=100)
    assert model.window_size == 100
    assert len(model.weights) == 101
    assert model.is_uniform
    assert model.max_weight == 1.0

def test_gc_bias_model_invalid_init():
    with pytest.raises(ValueError):
        GCBiasModel([1.0] * 100, 100)
    with pytest.raises(ValueError):
        GCBiasModel([1.0] * 101, 0)
    with pytest.raises(ValueError):
        GCBiasModel([-1.0] * 101, 100)

def test_gc_bias_model_get_weight():
    weights = [0.0] * 101
    weights[50] = 1.0
    model = GCBiasModel(weights, window_size=100)
    
    assert model.get_weight(0.5) == 1.0
    assert model.get_weight(0.496) == 1.0 # rounds to 0.50
    assert model.get_weight(0.0) == 0.0
    assert model.get_weight(1.0) == 0.0

def test_gc_bias_model_get_weight_for_sequence():
    weights = [0.0] * 101
    weights[50] = 1.0
    model = GCBiasModel(weights, window_size=100)
    
    # 50% GC
    assert model.get_weight_for_sequence("GCAT") == 1.0
    # 100% GC
    assert model.get_weight_for_sequence("GGCC") == weights[100]
    # 0% GC
    assert model.get_weight_for_sequence("AATT") == weights[0]
    
    # Empty/N
    assert model.get_weight_for_sequence("NNNN") == 1.0

def test_gc_bias_model_save_load(tmp_path):
    weights = np.random.random(101).tolist()
    window_size = 150
    model = GCBiasModel(weights, window_size)
    
    path = tmp_path / "gc_model.pickle"
    model.save(path)
    
    loaded_model = GCBiasModel.from_file(path)
    assert loaded_model.window_size == window_size
    np.testing.assert_allclose(loaded_model.weights, weights)

def test_uniform_gc_model():
    model = get_uniform_gc_model(200)
    assert model.window_size == 200
    assert model.is_uniform
    assert model.get_weight(0.3) == 1.0
