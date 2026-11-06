"""Frozen ecological arithmetic operating on precomputed response lookups.
The full archived lake replay also needs forcing, initialization and transport contracts.
"""
from .config import V16Config
from .production import project_two_species_conservative_production
from .state import simulate_grid_coupled_generator
