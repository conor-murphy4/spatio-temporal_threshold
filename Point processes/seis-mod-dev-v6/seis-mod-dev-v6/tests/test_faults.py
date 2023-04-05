import os
import unittest

import numpy as np
import pymc3 as pm

from seismod.geomodels import FaultModel


class TestFaultModelInitialization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path_faults = os.path.join(os.path.dirname(__file__), 'input_reference', 'fault_model.csv')

    def test_when_initializing_with_irregular_grid_size_then_fail(self):
        with self.assertRaisesRegex(ValueError, 'Grid spacing is not regular'):
            FaultModel(locations=[], throws=[], thickness=[], xnodes=[0, 1.5, 2], ynodes=[0, 1])

    def test_when_initializing_then_ratios_are_stored_as_throw_to_thickness(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        locs = [[1, 1]]
        throws, thick = [1], [2]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes)

        expected = np.zeros((1, 3, 3))
        expected[0, 1, 1] = .5
        np.testing.assert_equal(faults.ratio, expected)

    def test_when_initializing_from_csv_then_thickness_is_average_of_all_thickness_values(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        faults = FaultModel.from_csv(self.path_faults, xnodes=xnodes, ynodes=ynodes, dmin=0)

        expected = np.zeros((2, 3, 3))
        expected[:, 0, 0] = [1, .5]
        expected[0, 2, 2] = .5
        np.testing.assert_equal(faults.thickness, expected)

    def test_when_initializing_from_csv_then_throw_is_absolute_value_of_top_values_difference(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        faults = FaultModel.from_csv(self.path_faults, xnodes=xnodes, ynodes=ynodes, dmin=0)

        expected = np.zeros((2, 3, 3))
        expected[:, 0, 0] = [1, 2]
        expected[0, 2, 2] = 1
        np.testing.assert_equal(faults.thickness * faults.ratio, expected)

    def test_when_faults_are_closer_than_dmin_then_keep_only_first_fault(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        locs = [[1, 1], [1, 1.1], [2, 2]]
        throws, thick = [1, 1, 1], [1, 2, 3]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes, dmin=0.2)

        expected = np.zeros((1, 3, 3))
        expected[0, 1, 1] = 1
        expected[0, 2, 2] = 3
        np.testing.assert_equal(faults.thickness, expected)

    def test_when_multiple_locations_fall_in_same_grid_node_then_thickness_has_multiple_values_in_that_node(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        locs = [[1.1, 1.9], [0.9, 2.1]]
        throws, thick = [1, 1], [1, 2]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes, dmin=0)

        expected = np.zeros((2, 3, 3))
        expected[:, 1, 2] = [1, 2]
        np.testing.assert_equal(faults.thickness, expected)


class TestFaultModelTopographicGradientGridCreation(unittest.TestCase):
    def test_when_generating_topographic_gradient_grid_without_smoothing_then_throw_values_are_used(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        locs = [[1, 1]]
        throws, thick = [1], [3]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes)

        grid = faults.generate_topographic_gradients()
        expected = np.zeros((3, 3))
        expected[1, 1] = .5  # Throw divided by grid size (gradient approximation)
        np.testing.assert_equal(grid.attributes, expected)

    def test_when_generating_topographic_gradient_grid_without_max_ratio_then_max_throw_in_node_is_used(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        locs = [[1, 1], [1.1, 1], [1, 1.1]]
        throws, thick = [1, 3, 2], [3, 2, 4]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes, dmin=0)

        grid = faults.generate_topographic_gradients()
        expected = np.zeros((3, 3))
        expected[1, 1] = 1.5  # Throw divided by grid size (gradient approximation)
        np.testing.assert_equal(grid.attributes, expected)

    def test_when_generating_topographic_gradient_grid_with_max_ratio_then_max_throw_in_node_is_used(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        locs = [[1, 1], [1.1, 1], [1, 1.1]]
        throws, thick = [1, 3, 2], [3, 2, 4]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes, dmin=0)

        grid = faults.generate_topographic_gradients(rmax=1)
        expected = np.zeros((3, 3))
        expected[1, 1] = 1  # Throw divided by grid size (gradient approximation)
        np.testing.assert_equal(grid.attributes, expected)

    def test_when_generating_tensor_topographic_gradient_grid_with_max_ratio_then_max_throw_in_node_is_used(self):
        xnodes, ynodes = [0, 1, 2], [0, 1, 2]
        locs = [[1, 1], [1.1, 1], [1, 1.1]]
        throws, thick = [1, 3, 2], [3, 2, 4]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes, dmin=0)
        with pm.Model():
            rmax = pm.Uniform(name='a', testval=1)

        grid = faults.generate_topographic_gradients(rmax=rmax)
        expected = np.zeros((3, 3))
        expected[1, 1] = 1  # Throw divided by grid size (gradient approximation)
        np.testing.assert_equal(grid.attributes.tag.test_value, expected)

    def test_when_generating_topographic_gradient_grid_with_point_source_smoothing_then_grid_is_gaussian(self):
        xnodes, ynodes = np.arange(9), np.arange(9)
        locs = [[4, 4]]
        throws, thick = [2], [3]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes)

        grid = faults.generate_topographic_gradients(sigma=1.5)
        xcoord, ycoord = np.meshgrid(xnodes, ynodes, indexing='ij')
        expected = np.exp(-((xcoord - 4) ** 2 + (ycoord - 4) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        np.testing.assert_almost_equal(grid.attributes, expected)

    def test_when_generating_tensor_topographic_gradient_grid_with_point_source_smoothing_then_grid_is_gaussian(self):
        xnodes, ynodes = np.arange(9), np.arange(9)
        locs = [[4, 4]]
        throws, thick = [2], [3]
        faults = FaultModel(locations=locs, throws=throws, thickness=thick, xnodes=xnodes, ynodes=ynodes)

        with pm.Model():
            sigma = pm.Uniform(name='a', upper=2, testval=1.5)
        grid = faults.generate_topographic_gradients(sigma=sigma)
        xcoord, ycoord = np.meshgrid(xnodes, ynodes, indexing='ij')
        expected = np.exp(-((xcoord - 4) ** 2 + (ycoord - 4) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        np.testing.assert_almost_equal(grid.attributes.tag.test_value, expected)


if __name__ == '__main__':
    unittest.main()
