import unittest
from datetime import datetime

import numpy as np

from seismod.geomodels import EarthquakeCatalogue


class TestEarthquakeCatalogueInitialization(unittest.TestCase):
    def test_when_initializing_without_parameters_then_create_empty_catalogue(self):
        catalogue = EarthquakeCatalogue()
        self.assertTupleEqual(catalogue.nevents, (0,))

    def test_when_initializing_with_magnitude_only_then_magnitude_has_same_values(self):
        catalogue = EarthquakeCatalogue(magnitude=[1., 2.])
        np.testing.assert_equal(catalogue.magnitude, [1., 2.])

    def test_when_initializing_with_magnitude_only_then_location_has_default_values(self):
        catalogue = EarthquakeCatalogue(magnitude=[1., 2.])
        np.testing.assert_equal(catalogue.location, [[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])

    def test_when_initializing_with_magnitude_only_then_time_has_default_values(self):
        catalogue = EarthquakeCatalogue(magnitude=[1., 2.])
        np.testing.assert_equal(catalogue.time, np.full(2, 'NaT', dtype='M8[s]'))

    def test_when_initializing_with_magnitude_only_then_trial_has_default_values(self):
        catalogue = EarthquakeCatalogue(magnitude=[1., 2.])
        np.testing.assert_equal(catalogue.trial, [0, 0])

    def test_when_initializing_with_2d_catalogue_then_trial_tracks_second_dimension(self):
        catalogue = EarthquakeCatalogue(magnitude=[[1., 2.], [2., 1.]])
        np.testing.assert_equal(catalogue.trial, [[0, 1], [0, 1]])

    def test_when_initializing_with_2d_catalogue_and_explicit_trial_then_fail(self):
        with self.assertRaisesRegex(ValueError, '2D-catalogue with trial'):
            EarthquakeCatalogue(trial=[[0, 1], [2, 3]])

    def test_when_initializing_with_more_than_2d_catalogue_then_fail(self):
        with self.assertRaisesRegex(ValueError, 'Too many dimensions'):
            EarthquakeCatalogue(magnitude=[[[1.], [2.]], [[2.], [1.]]])

    def test_when_initializing_with_unknown_keywords_then_fail(self):
        with self.assertRaisesRegex(KeyError, 'Unrecognized arguments'):
            EarthquakeCatalogue(unknown=[0, 2])

    def test_when_initializing_with_mismatched_argument_shapes_then_fail(self):
        with self.assertRaisesRegex(ValueError, 'Mismatched dimensions'):
            EarthquakeCatalogue(magnitude=[0, 2], trial=[0, 1, 2])

    def test_when_initializing_with_location_with_more_than_3_coordinates_then_fail(self):
        with self.assertRaisesRegex(ValueError, 'with dimension 3'):
            EarthquakeCatalogue(location=[[0, 0, 0, 0], [0, 0, 0, 0]])


class TestEarthquakeCatalogueProperties(unittest.TestCase):
    def test_when_computing_inter_event_distance_then_result_is_cartesian_distance(self):
        catalogue = EarthquakeCatalogue(location=[[0, 0, 0], [1, 1, 1]])
        self.assertAlmostEqual(catalogue.dr[1, 0], np.sqrt(3))

    def test_when_computing_inter_event_distance_for_2d_catalogue_then_result_is_computed_for_each_catalogue(self):
        catalogue = EarthquakeCatalogue(location=[[[0, 0, 0], [2, 2, 2]], [[1, 1, 1], [4, 4, 4]]])
        np.testing.assert_allclose(catalogue.dr[1, 0], [np.sqrt(3), np.sqrt(12)])

    def test_when_computing_inter_event_time_then_result_is_in_days_wrt_first_dimension(self):
        catalogue = EarthquakeCatalogue(time=[datetime(2020, 1, 1), datetime(2020, 2, 1)])
        self.assertAlmostEqual(catalogue.dt[1, 0], 31)

    def test_when_computing_inter_event_time_for_2d_catalogue_then_result_is_computed_for_each_catalogue(self):
        catalogue = EarthquakeCatalogue(time=[[datetime(2020, 1, 1), datetime(2020, 2, 1)],
                                              [datetime(2020, 2, 2), datetime(2020, 3, 2)]])
        np.testing.assert_allclose(catalogue.dt[1, 0], [32, 30])

    def test_when_computing_inter_event_distance_then_matrix_is_symmetric(self):
        catalogue = EarthquakeCatalogue(location=np.random.random((20, 3)))
        np.testing.assert_allclose(catalogue.dr, catalogue.dr.T)

    def test_when_computing_inter_event_time_then_matrix_is_antisymmetric(self):
        delta_in_seconds = np.random.randint(365 * 2 * 24 * 3600, size=20)
        catalogue = EarthquakeCatalogue(time=np.datetime64('2020-01-01 00:00:00') + delta_in_seconds)
        np.testing.assert_allclose(catalogue.dt, -catalogue.dt.T)

    def test_when_computing_inter_event_distance_then_diagonal_is_zero(self):
        catalogue = EarthquakeCatalogue(location=np.random.random((20, 3)))
        np.testing.assert_allclose(catalogue.dr.diagonal(), 0)

    def test_when_computing_inter_event_distance_for_2d_catalogue_then_each_catalogue_submatrix_is_symmetric(self):
        catalogue = EarthquakeCatalogue(location=np.random.random((20, 2, 3)))
        np.testing.assert_allclose(catalogue.dr[..., 0], catalogue.dr[..., 0].T)
        np.testing.assert_allclose(catalogue.dr[..., 1], catalogue.dr[..., 1].T)

    def test_when_computing_inter_event_time_for_2d_catalogue_then_each_catalogue_submatrix_is_antisymmetric(self):
        delta_in_seconds = np.random.randint(365 * 2 * 24 * 3600, size=(20, 2))
        catalogue = EarthquakeCatalogue(time=np.datetime64('2020-01-01 00:00:00') + delta_in_seconds)
        np.testing.assert_allclose(catalogue.dt[..., 0], -catalogue.dt[..., 0].T)
        np.testing.assert_allclose(catalogue.dt[..., 1], -catalogue.dt[..., 1].T)

    def test_when_computing_inter_event_distance_for_2d_catalogue_then_diagonal_is_zero(self):
        catalogue = EarthquakeCatalogue(location=np.random.random((20, 2, 3)))
        np.testing.assert_allclose(catalogue.dr.diagonal(), 0)

    def test_when_getting_Mmax_in_empty_catalogue_then_Mmax_is_NaN(self):
        catalogue = EarthquakeCatalogue()
        self.assertTrue(np.isnan(catalogue.Mmax))

    def test_when_getting_Mmax_in_1d_catalogue_then_maximum_of_whole_catalogue_is_returned(self):
        catalogue = EarthquakeCatalogue(magnitude=[2, 3, 4.5, 4])
        self.assertEqual(catalogue.Mmax, 4.5)

    def test_when_getting_Mmax_in_2d_catalogue_then_maximum_per_catalogue_is_returned(self):
        catalogue = EarthquakeCatalogue(magnitude=[[2, 4], [4.5, 3]])
        np.testing.assert_equal(catalogue.Mmax, [4.5, 4])


class TestEarthquakeCatalogueSortingAndSelecting(unittest.TestCase):
    def test_when_sorting_by_time_then_magnitudes_are_sorted_accordingly(self):
        catalogue = EarthquakeCatalogue(magnitude=[2, 3], time=[datetime(2020, 1, 2), datetime(2020, 1, 1)])
        catalogue.sort()
        np.testing.assert_equal(catalogue.magnitude, [3, 2])

    def test_when_sorting_2d_catalogue_by_time_then_magnitudes_in_each_catalogue_are_sorted_accordingly(self):
        catalogue = EarthquakeCatalogue(magnitude=[[2, 3], [4, 1]], time=[[datetime(2020, 1, 2), datetime(2020, 1, 3)],
                                                                          [datetime(2020, 1, 1), datetime(2020, 1, 4)]])
        catalogue.sort()
        np.testing.assert_equal(catalogue.magnitude, [[4, 3], [2, 1]])

    def test_when_selecting_by_trial_in_1d_catalogue_then_selected_trials_are_kept_in_ascending_order(self):
        catalogue = EarthquakeCatalogue(magnitude=[1, 2, 3, 4], trial=[0, 1, 2, 3])
        catalogue.select(trial=[0, 3, 2])
        np.testing.assert_equal(catalogue.magnitude, [1, 3, 4])

    def test_when_selecting_by_trial_with_non_existing_trial_then_select_existing_only_and_warn(self):
        catalogue = EarthquakeCatalogue(magnitude=[1, 2, 3, 4], trial=[0, 1, 2, 3])
        with self.assertLogs('geomodels', 'WARNING'):
            catalogue.select(trial=[0, 3, 2, 4])
        np.testing.assert_equal(catalogue.magnitude, [1, 3, 4])

    def test_when_selecting_by_trial_in_2d_catalogue_then_selected_trials_are_kept_in_ascending_order(self):
        catalogue = EarthquakeCatalogue(magnitude=[[1, 2, 3, 4]])
        catalogue.select(trial=[0, 3, 2])
        np.testing.assert_equal(catalogue.magnitude, [[1, 3, 4]])

    def test_when_selecting_then_minimum_is_inclusive_and_maximum_is_exclusive(self):
        catalogue = EarthquakeCatalogue(magnitude=[1, 2, 3, 4])
        catalogue.select(Mmin=2, Mmax=4)
        np.testing.assert_equal(catalogue.magnitude, [2, 3])

    def test_when_selecting_with_mask_in_2d_catalogue_then_warn_and_catalogue_is_flattened(self):
        catalogue = EarthquakeCatalogue(magnitude=[[1, 2], [3, 4]])
        with self.assertLogs('geomodels', 'WARNING'):
            catalogue.select(Mmin=2, Mmax=4)
        np.testing.assert_equal(catalogue.magnitude, [2, 3])
        np.testing.assert_equal(catalogue.trial, [1, 0])


if __name__ == '__main__':
    unittest.main()
