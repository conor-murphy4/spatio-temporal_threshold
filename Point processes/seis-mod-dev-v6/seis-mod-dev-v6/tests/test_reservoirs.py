import os
import unittest
import datetime
from unittest import mock

import numpy as np
import pymc3 as pm
from scipy import interpolate

from seismod.geomodels import ReservoirGrid, ReservoirGridEpochs


class TestReservoirGridInitialization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path_reservoir = os.path.join(os.path.dirname(__file__), 'input_reference', 'reservoir.csv')

    def test_when_initializing_with_mismatched_argument_shapes_then_fail(self):
        with self.assertRaisesRegex(ValueError, 'Mismatched dimensions'):
            ReservoirGrid(names=['a'], xnodes=[0, 1], ynodes=[0, 1, 2], attributes=np.zeros((1, 2, 2)))

    def test_when_initializing_with_irregular_grid_size_then_fail(self):
        with self.assertRaisesRegex(ValueError, 'Grid spacing is not regular'):
            ReservoirGrid(names=['a'], xnodes=[0, 1.5, 2], ynodes=[0, 1], attributes=np.zeros((1, 3, 2)))

    def test_when_initializing_with_one_node_then_fail(self):
        with self.assertRaisesRegex(IndexError, 'index 0 is out of bounds'):
            ReservoirGrid(names=['a'], xnodes=[0], ynodes=[0], attributes=np.zeros((1, 1, 1)))

    def test_when_initializing_from_csv_then_nodes_are_read_correctly_and_missing_node_is_filled_in(self):
        grid = ReservoirGrid.from_csv(self.path_reservoir, fill=-1)
        np.testing.assert_equal(grid.attributes, [[1, 1, 1], [1, -1, 1], [1, 1, 1]])

    def test_when_initializing_from_csv_then_missing_mask_matches_missing_node(self):
        grid = ReservoirGrid.from_csv(self.path_reservoir, fill=-1)
        np.testing.assert_equal(grid.mask, [[[False, False, False], [False, True, False], [False, False, False]]])

    def test_when_initializing_from_csv_then_attribute_names_are_read_correctly(self):
        grid = ReservoirGrid.from_csv(self.path_reservoir, fill=-1)
        self.assertTupleEqual(grid.names, ('attribute',))

    def test_when_initializing_then_grid_spacing_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a'], xnodes=[0, 1.5, 3], ynodes=[0, 1.5], attributes=np.zeros((1, 3, 2)))
        self.assertAlmostEqual(grid.dx, 1.5)
        self.assertAlmostEqual(grid.dy, 1.5)

    def test_when_initializing_then_minimum_easting_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a'], xnodes=[0, 1.5, 3], ynodes=[2, 3.5], attributes=np.zeros((1, 3, 2)))
        self.assertEqual(grid.minx, 0)

    def test_when_initializing_then_maximum_easting_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a'], xnodes=[0, 1.5, 3], ynodes=[2, 3.5], attributes=np.zeros((1, 3, 2)))
        self.assertEqual(grid.maxx, 3)

    def test_when_initializing_then_minimum_northing_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a'], ynodes=[0, 1.5, 3], xnodes=[2, 3.5], attributes=np.zeros((1, 2, 3)))
        self.assertEqual(grid.miny, 0)

    def test_when_initializing_then_maximum_northing_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a'], ynodes=[0, 1.5, 3], xnodes=[2, 3.5], attributes=np.zeros((1, 2, 3)))
        self.assertEqual(grid.maxy, 3)

    def test_when_initializing_then_number_of_easting_nodes_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a'], xnodes=[0, 1, 2], ynodes=[5, 6], attributes=np.zeros((1, 3, 2)))
        self.assertEqual(grid.nx, 3)

    def test_when_initializing_then_number_of_northing_nodes_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a'], ynodes=[0, 1, 2], xnodes=[5, 6], attributes=np.zeros((1, 2, 3)))
        self.assertEqual(grid.ny, 3)

    def test_when_initializing_then_number_of_attributes_is_calculated_correctly(self):
        grid = ReservoirGrid(names=['a', 'b'], ynodes=[0, 1], xnodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        self.assertEqual(grid.nattrs, 2)

    def test_when_initializing_then_easting_nodes_is_unmodified(self):
        grid = ReservoirGrid(names=['a'], xnodes=[0, 1, 2], ynodes=[5, 6], attributes=np.zeros((1, 3, 2)))
        np.testing.assert_equal(grid.xnodes, [0, 1, 2])

    def test_when_initializing_then_northing_nodes_is_unmodified(self):
        grid = ReservoirGrid(names=['a'], ynodes=[0, 1, 2], xnodes=[5, 6], attributes=np.zeros((1, 2, 3)))
        np.testing.assert_equal(grid.ynodes, [0, 1, 2])

    def test_when_initializing_then_attribute_names_are_unmodified(self):
        grid = ReservoirGrid(names=['a', 'b'], ynodes=[0, 1], xnodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        self.assertTupleEqual(grid.names, ('a', 'b'))

    def test_when_initializing_with_non_nan_fill_value_then_filled_value_mask_matches_fill_value(self):
        grid = ReservoirGrid(names=['a'], ynodes=[0, 1], xnodes=[0, 1], attributes=[[[0, 1], [1, 0]]], fill=1)
        np.testing.assert_equal(grid.mask, [[[False, True], [True, False]]])

    def test_when_initializing_with_nan_fill_value_then_filled_value_mask_matches_fill_value(self):
        nan = np.nan
        grid = ReservoirGrid(names=['a'], ynodes=[0, 1], xnodes=[0, 1], attributes=[[[nan, 1], [1, nan]]], fill=nan)
        np.testing.assert_equal(grid.mask, [[[True, False], [False, True]]])


class TestReservoirGridFilters(unittest.TestCase):
    def test_when_snapping_then_fill_in_values_are_replaced_with_spatial_nearest_neighbors(self):
        xnodes = [0, 1, 2]
        ynodes = [0, 1, 2]
        attributes = [[[0, 1, np.nan], [2, 3, 4], [np.nan, 5, 6]]]
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=np.nan)
        grid.snap()
        np.testing.assert_equal(grid['a'], [[0, 1, 1], [2, 3, 4], [2, 5, 6]])

    def test_when_gaussian_filtering_point_source_then_result_is_gaussian(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5, 6]
        attributes = np.zeros((1, 7, 7))
        attributes[0, 3, 3] = 1
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=np.nan)
        grid.gaussian_filter(1.5)

        xcoord, ycoord = np.meshgrid(xnodes, ynodes, indexing='ij')
        expected = np.exp(-((xcoord - 3) ** 2 + (ycoord - 3) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        np.testing.assert_allclose(grid.attributes, expected)

    def test_when_gaussian_filtering_off_center_point_source_then_result_is_off_center_gaussian(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5]
        attributes = np.zeros((1, 7, 6))
        attributes[0, 3, 3] = 1
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=np.nan)
        grid.gaussian_filter(1.5)

        xcoord, ycoord = np.meshgrid(np.arange(7), np.arange(7), indexing='ij')
        expected = np.exp(-((xcoord - 3) ** 2 + (ycoord - 3) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        # Symmetrical boundary conditions
        expected[:, -2] += expected[:, -1]
        # Slicing due to how padding works in even sized dimensions
        np.testing.assert_allclose(grid.attributes, expected[:, :-1])

    def test_when_gaussian_filtering_tensor_point_source_then_result_is_gaussian_tensor(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5, 6]
        attributes = np.zeros((1, 7, 7))
        attributes[0, 3, 3] = 1
        with pm.Model():
            attr = pm.Uniform(name='a', shape=(1, 7, 7), testval=attributes)
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attr, fill=np.nan)
        grid.gaussian_filter(1.5)

        xcoord, ycoord = np.meshgrid(xnodes, ynodes, indexing='ij')
        expected = np.exp(-((xcoord - 3) ** 2 + (ycoord - 3) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        np.testing.assert_allclose(grid.attributes.tag.test_value, expected)

    def test_when_custom_gaussian_filtering_point_source_then_result_is_gaussian(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5, 6]
        attributes = np.zeros((1, 7, 7))
        attributes[0, 3, 3] = 1
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=np.nan)
        grid.adapted_gaussian_filter(1.5)

        xcoord, ycoord = np.meshgrid(xnodes, ynodes, indexing='ij')
        expected = np.exp(-((xcoord - 3) ** 2 + (ycoord - 3) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        np.testing.assert_allclose(grid.attributes, expected)

    def test_when_custom_gaussian_filtering_off_center_point_source_then_result_is_off_center_gaussian(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5]
        attributes = np.zeros((1, 7, 6))
        attributes[0, 3, 3] = 1
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=np.nan)
        grid.adapted_gaussian_filter(1.5)

        xcoord, ycoord = np.meshgrid(np.arange(7), np.arange(7), indexing='ij')
        expected = np.exp(-((xcoord - 3) ** 2 + (ycoord - 3) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        # Additional factor due to kernel overlap (because of custom filtering method)
        expected[3, 3] += np.sum(expected[:, -1])
        # Slicing due to how padding works in even sized dimensions
        np.testing.assert_allclose(grid.attributes, expected[:, :-1])

    def test_when_custom_gaussian_filtering_tensor_point_source_then_result_is_gaussian_tensor(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5, 6]
        attributes = np.zeros((1, 7, 7))
        attributes[0, 3, 3] = 1
        with pm.Model():
            attr = pm.Uniform(name='a', shape=(1, 7, 7), testval=attributes)
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attr, fill=np.nan)
        grid.adapted_gaussian_filter(1.5)

        xcoord, ycoord = np.meshgrid(xnodes, ynodes, indexing='ij')
        expected = np.exp(-((xcoord - 3) ** 2 + (ycoord - 3) ** 2) / (2 * 1.5 ** 2))
        expected /= np.sum(expected)
        np.testing.assert_allclose(grid.attributes.tag.test_value, expected)

    def test_when_median_filtering_with_even_kernel_size_then_fail(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5, 6]
        attributes = np.zeros((1, 7, 7))
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=np.nan)
        with self.assertRaisesRegex(ValueError, 'must be odd'):
            grid.median_filter(2)

    def test_when_median_filtering_with_negative_kernel_size_then_fail(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5, 6]
        attributes = np.zeros((1, 7, 7))
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=np.nan)
        with self.assertRaisesRegex(ValueError, 'must be positive'):
            grid.median_filter(-2)

    def test_when_median_filtering_then_ignore_masked_values(self):
        xnodes = [0, 1, 2, 3, 4, 5, 6]
        ynodes = [0, 1, 2, 3, 4, 5, 6]
        attributes = np.zeros((1, 7, 7))
        attributes[:, 2:-2, 2:-2] = 1
        attributes[:, 3, 3] = 1
        grid = ReservoirGrid(names=['a'], ynodes=ynodes, xnodes=xnodes, attributes=attributes, fill=0)
        grid.median_filter(3)

        expected = np.zeros((7, 7))
        expected[2:-2, 2:-2] = 1
        np.testing.assert_allclose(grid.attributes, expected)


class TestReservoirGridDataAccess(unittest.TestCase):
    def test_when_accessing_specific_attribute_data_then_output_is_corresponding_index_attribute(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGrid(names=['a', 'b'], ynodes=[0, 1], xnodes=[0, 1], attributes=attributes)
        np.testing.assert_allclose(grid['b'], [[4, 5], [6, 7]])

    def test_when_accessing_location_data_with_1d_array_then_output_is_squeezed(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGrid(names=['a', 'b'], ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        np.testing.assert_allclose(grid.value_at_x(np.array([0, 1.5])), [1, 5])

    def test_when_accessing_exact_location_data_then_output_matches_attributes_exactly(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGrid(names=['a', 'b'], ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        np.testing.assert_allclose(grid.value_at_x(np.array([[0, 1.5], [0, 1.5], [0, 0]])), [[1, 1, 0], [5, 5, 4]])

    def test_when_accessing_location_data_outside_grid_then_output_correspond_to_location_clipped_to_grid(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGrid(names=['a', 'b'], ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        np.testing.assert_allclose(grid.value_at_x(np.array([[-3, -3], [4, -3]])), [[0, 2], [4, 6]])

    def test_when_accessing_non_exact_location_data_then_output_matches_attributes_with_closest_coordinates(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGrid(names=['a', 'b'], ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        np.testing.assert_allclose(grid.value_at_x(np.array([[.75, .76], [.76, .76]])), [[1, 3], [5, 7]])


class TestReservoirGridEpochsInitialization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path_reservoir = os.path.join(os.path.dirname(__file__), 'input_reference', 'reservoir_epochs.csv')

    def test_when_initializing_with_non_constant_epoch_spacing_then_fail(self):
        epochs = ['P2010_01_01', 'P2011_01_01', 'P2012_01_02']
        with self.assertRaisesRegex(ValueError, 'Non-constant epoch spacing'):
            ReservoirGridEpochs(names=epochs, xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((3, 2, 2)))

    def test_when_initializing_with_non_annual_epoch_spacing_then_fail(self):
        epochs = ['P2010_01_01', 'P2010_02_01', 'P2010_03_01']
        with self.assertRaisesRegex(NotImplementedError, 'Only annual epochs are currently supported'):
            ReservoirGridEpochs(names=epochs, xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((3, 2, 2)))

    def test_when_initializing_with_non_constant_prefix_then_fail(self):
        epochs = ['P2010_01_01', 'p2011_01_01']
        with self.assertRaisesRegex(ValueError, 'Multiple formats for epoch labels'):
            ReservoirGridEpochs(names=epochs, xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((2, 2, 2)))

    def test_when_initializing_with_multiple_date_formats_then_fail(self):
        epochs = ['P2010_01_01', 'P2011']
        with self.assertRaisesRegex(ValueError, 'Multiple formats for epoch labels'):
            ReservoirGridEpochs(names=epochs, xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((2, 2, 2)))

    def test_when_initializing_from_csv_then_nodes_are_read_correctly_and_missing_node_is_filled_in(self):
        grid = ReservoirGridEpochs.from_csv(self.path_reservoir, fill=-1)
        np.testing.assert_equal(grid.attributes, [[[1, 1, 1], [1, -1, 1], [1, 1, 1]],
                                                  [[2, 2, 2], [2, -1, 2], [2, 2, 2]]])

    def test_when_initializing_from_csv_then_missing_mask_matches_missing_node(self):
        grid = ReservoirGridEpochs.from_csv(self.path_reservoir, fill=-1)
        np.testing.assert_equal(grid.mask, [[[False, False, False], [False, True, False], [False, False, False]],
                                            [[False, False, False], [False, True, False], [False, False, False]]])

    def test_when_initializing_from_csv_then_epochs_are_read_correctly(self):
        grid = ReservoirGridEpochs.from_csv(self.path_reservoir, fill=-1)
        np.testing.assert_equal(grid.epochs, np.array(['2020-10-01', '2021-10-01'], dtype='M8[s]'))

    def test_when_using_only_years_in_label_then_all_epochs_are_set_on_january_first(self):
        labels = ['P2010', 'P2011']
        grid = ReservoirGridEpochs(names=labels, xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        np.testing.assert_equal(grid.epochs, np.array(['2010-01-01', '2011-01-01'], dtype='M8[D]'))

    def test_when_using_year_month_date_in_label_then_all_epochs_match_day_and_month(self):
        labels = ['P2010_09_13', 'P2011_09_13']
        grid = ReservoirGridEpochs(names=labels, xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        np.testing.assert_equal(grid.epochs, np.array(['2010-09-13', '2011-09-13'], dtype='M8[D]'))

    def test_when_using_arbitrary_prefix_in_label_then_prefix_is_handled_correctly(self):
        labels = [';a[strange&prefix2010_01_01', ';a[strange&prefix2011_01_01']
        grid = ReservoirGridEpochs(names=labels, xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        np.testing.assert_equal(grid.epochs, np.array(['2010-01-01', '2011-01-01'], dtype='M8[D]'))


class TestReservoirGridEpochsDataAccess(unittest.TestCase):
    def test_when_accessing_attribute_by_epoch_as_datetime_then_output_is_corresponding_attribute_data(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGridEpochs(names=['2020', '2021'], ynodes=[0, 1], xnodes=[0, 1], attributes=attributes)
        np.testing.assert_allclose(grid[datetime.datetime(2021, 1, 1)], [[4, 5], [6, 7]])

    def test_when_accessing_attribute_by_epoch_as_numpy_datetime_then_output_is_corresponding_attribute_data(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGridEpochs(names=['2020', '2021'], ynodes=[0, 1], xnodes=[0, 1], attributes=attributes)
        np.testing.assert_allclose(grid[np.datetime64('2021-01-01')], [[4, 5], [6, 7]])

    def test_when_accessing_attribute_by_attribute_name_then_output_is_corresponding_attribute_data(self):
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGridEpochs(names=['dV2020', 'dV2021'], ynodes=[0, 1], xnodes=[0, 1], attributes=attributes)
        np.testing.assert_allclose(grid['dV2021'], [[4, 5], [6, 7]])

    def test_when_accessing_epoch_subset_then_start_and_end_date_are_included(self):
        labels = ['2020', '2021', '2022', '2023']
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1], xnodes=[0, 1], attributes=np.zeros((4, 2, 2)))
        subset = grid.period_epochs(datetime.datetime(2021, 1, 1), datetime.datetime(2022, 1, 1))
        np.testing.assert_equal(subset, np.array(['2021-01-01', '2022-01-01'], dtype='M8[D]'))

    def test_when_accessing_epoch_subset_with_start_outside_epochs_then_fail(self):
        labels = ['2020', '2021']
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1], xnodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        with self.assertRaisesRegex(ValueError, 'Period start or end not present'):
            grid.period_epochs(datetime.datetime(2019, 1, 1), datetime.datetime(2021, 1, 1))

    def test_when_accessing_epoch_subset_with_end_outside_epochs_then_fail(self):
        labels = ['2020', '2021']
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1], xnodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        with self.assertRaisesRegex(ValueError, 'Period start or end not present'):
            grid.period_epochs(datetime.datetime(2020, 1, 1), datetime.datetime(2022, 1, 1))

    def test_when_accessing_attribute_subset_then_start_and_end_date_attribute_data_are_included(self):
        labels = ['2020', '2021', '2022', '2023']
        attributes = np.arange(4 * 2 * 2).reshape((4, 2, 2))
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1], xnodes=[0, 1], attributes=attributes)
        subattributes = grid.period_attributes(datetime.datetime(2021, 1, 1), datetime.datetime(2022, 1, 1))
        np.testing.assert_equal(subattributes, np.arange(2 * 2, 3 * 2 * 2).reshape((2, 2, 2)))

    def test_when_accessing_attribute_subset_with_start_outside_epochs_then_fail(self):
        labels = ['2020', '2021']
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1], xnodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        with self.assertRaisesRegex(ValueError, 'Period start or end not present'):
            grid.period_attributes(datetime.datetime(2019, 1, 1), datetime.datetime(2021, 1, 1))

    def test_when_accessing_attribute_subset_with_end_outside_epochs_then_fail(self):
        labels = ['2020', '2021']
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1], xnodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        with self.assertRaisesRegex(ValueError, 'Period start or end not present'):
            grid.period_attributes(datetime.datetime(2020, 1, 1), datetime.datetime(2022, 1, 1))

    def test_when_accessing_attribute_at_exact_time_and_location_then_output_matches_attributes_exactly(self):
        labels = ['2020', '2021']
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        time, loc = np.array(['2021-01-01', '2020-01-01'], dtype='M8[D]'), np.array([[0, 1.5], [1.5, 1.5]])
        np.testing.assert_equal(grid.value_at_tx(time, loc), [5, 3])

    def test_when_accessing_attribute_at_arbitrary_time_with_2_epochs_then_output_is_linearly_interpolated(self):
        labels = ['2020', '2021']
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        time, loc = np.array(['2020-01-30T06:30:00'], dtype='M8[s]'), np.array([[0, 0]])
        np.testing.assert_allclose(grid.value_at_tx(time, loc), [((29 + 6.5 / 24) / 366) * 4])

    def test_when_accessing_attribute_at_time_and_location_with_scalar_input_then_output_is_scalar(self):
        labels = ['2020', '2021']
        attributes = [[[0, 1], [2, 3]], [[4, 5], [6, 7]]]
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        time, loc = np.datetime64('2020-01-30T06:30:00'), np.array([0, 0])
        self.assertAlmostEqual(grid.value_at_tx(time, loc), ((29 + 6.5 / 24) / 366) * 4)

    def test_when_accessing_attribute_at_arbitrary_time_then_output_is_cubic_spline_interpolated(self):
        labels = ['2020', '2021', '2022', '2023']
        attributes = np.zeros((4, 2, 2))
        attributes[:, 0, 0] = 4 * np.arange(4) ** 3 - 2 * np.arange(4) ** 2 + np.arange(4)
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        time, loc = np.array(['2021-03-01T06:30:00'], dtype='M8[s]'), np.array([[0, 0]])
        dt = 1 + (59 + 6.5 / 24) / 365
        np.testing.assert_allclose(grid.value_at_tx(time, loc), [4 * dt ** 3 - 2 * dt ** 2 + dt], rtol=1e-4)

    def test_when_accessing_attribute_at_time_and_location_with_too_large_array_then_output_is_still_possible(self):
        labels = ['2020', '2021', '2022', '2023']
        attributes = np.zeros((4, 2, 2))
        attributes[:, 0, 0] = 4 * np.arange(4) ** 3 - 2 * np.arange(4) ** 2 + np.arange(4)
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        time, loc = np.array([['2021-03-01T06:30:00']], dtype='M8[s]'), np.array([[[0, 0]]])
        dt = 1 + (59 + 6.5 / 24) / 365

        result = [4 * dt ** 3 - 2 * dt ** 2 + dt]
        array_result = np.array([[result, result], [result, result]])
        with mock.patch.object(interpolate.CubicSpline, '__call__') as mock_interp:
            mock_interp.side_effect = [MemoryError(), array_result]
            np.testing.assert_allclose(grid.value_at_tx(time, loc), [result], rtol=1e-4)

    def test_when_accessing_attribute_at_time_and_location_with_too_large_array_and_single_dimension_then_fail(self):
        labels = ['2020', '2021', '2022', '2023']
        attributes = np.zeros((4, 2, 2))
        attributes[:, 0, 0] = 4 * np.arange(4) ** 3 - 2 * np.arange(4) ** 2 + np.arange(4)
        grid = ReservoirGridEpochs(names=labels, ynodes=[0, 1.5], xnodes=[0, 1.5], attributes=attributes)
        time, loc = np.array(['2021-03-01T06:30:00'], dtype='M8[s]'), np.array([[0, 0]])
        dt = 1 + (59 + 6.5 / 24) / 365

        result = [4 * dt ** 3 - 2 * dt ** 2 + dt]
        array_result = np.array([[result, result], [result, result]])
        with mock.patch.object(interpolate.CubicSpline, '__call__') as mock_interp:
            mock_interp.side_effect = [MemoryError(), array_result]
            with self.assertRaises(MemoryError):
                grid.value_at_tx(time, loc)


class TestReservoirGridNodeRemoval(unittest.TestCase):
    def test_when_all_inputs_are_none_then_run_without_modifying_anything(self):
        grid1, grid2 = None, None
        ReservoirGrid.match_nodes(grid1, grid2)
        self.assertIsNone(grid1)
        self.assertIsNone(grid2)

    def test_when_removing_nodes_then_keep_only_common_nodes(self):
        grid1 = ReservoirGrid(names=['a'], xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((1, 2, 2)))
        grid2 = ReservoirGrid(names=['b'], xnodes=[0, 1, 2], ynodes=[0, 1, 2], attributes=np.zeros((1, 3, 3)))

        ReservoirGrid.match_nodes(grid1, grid2)
        np.testing.assert_equal(grid2.xnodes, grid1.xnodes)
        np.testing.assert_equal(grid1.xnodes, [0, 1])
        np.testing.assert_equal(grid2.ynodes, grid1.ynodes)
        np.testing.assert_equal(grid1.ynodes, [0, 1])

    def test_when_removing_nodes_from_grids_with_different_node_size_then_fail(self):
        grid1 = ReservoirGrid(names=['a'], xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((1, 2, 2)))
        grid2 = ReservoirGrid(names=['b'], xnodes=[0, .5, 1], ynodes=[0, .5, 1], attributes=np.zeros((1, 3, 3)))

        with self.assertRaisesRegex(ValueError, 'grids with different grid spacing'):
            ReservoirGrid.match_nodes(grid1, grid2)

    def test_when_removing_nodes_then_do_not_modify_attribute_names(self):
        grid1 = ReservoirGrid(names=['a', 'b'], xnodes=[0, 1], ynodes=[0, 1], attributes=np.zeros((2, 2, 2)))
        grid2 = ReservoirGrid(names=['c', 'd'], xnodes=[0, 1, 2], ynodes=[0, 1, 2], attributes=np.zeros((2, 3, 3)))

        ReservoirGrid.match_nodes(grid1, grid2)
        self.assertTupleEqual(grid1.names, ('a', 'b'))
        self.assertTupleEqual(grid2.names, ('c', 'd'))

    def test_when_removing_nodes_then_selected_attribute_values_match_nodes(self):
        attributes1 = np.random.random((2, 2, 2))
        attributes2 = np.random.random((2, 3, 3))
        grid1 = ReservoirGrid(names=['a', 'b'], xnodes=[0, 1], ynodes=[0, 1], attributes=attributes1)
        grid2 = ReservoirGrid(names=['c', 'd'], xnodes=[0, 1, 2], ynodes=[0, 1, 2], attributes=attributes2)

        ReservoirGrid.match_nodes(grid1, grid2)
        np.testing.assert_equal(grid1.attributes, attributes1)
        np.testing.assert_equal(grid2.attributes, attributes2[:, :2, :2])


if __name__ == '__main__':
    unittest.main()
