import unittest
from unittest import mock

import numpy as np

from seismod.geomodels import ThinSheetModel, ReservoirGrid, ReservoirGridEpochs, FaultModel


class TestThinsheetInitialization(unittest.TestCase):
    def test_when_initializing_without_paths_then_initialize_empty_thinsheet(self):
        with mock.patch.object(ReservoirGrid, 'from_csv') as mock_rgrid:
            with mock.patch.object(ReservoirGridEpochs, 'from_csv') as mock_rgride:
                with mock.patch.object(FaultModel, 'from_csv') as mock_faults:
                    ThinSheetModel(paths=dict())
        mock_rgrid.assert_not_called()
        mock_rgride.assert_not_called()
        mock_faults.assert_not_called()

    def test_when_initializing_with_thickness_and_median_filter_then_initialize_filtered_thickness(self):
        with mock.patch.multiple(ReservoirGrid, from_csv=mock.DEFAULT, match_nodes=mock.DEFAULT) as mock_rgrid:
            thickness_grid = mock.MagicMock()
            mock_rgrid['from_csv'].return_value = thickness_grid
            ThinSheetModel(paths=dict(thickness='a', compressibility='b'), kernels=dict(thickness=3))
        self.assertEqual(mock_rgrid['from_csv'].call_count, 2)
        thickness_grid.median_filter.assert_called_once_with(3)
        mock_rgrid['match_nodes'].assert_called_once()

    def test_when_initializing_with_thickness_and_gaussian_filter_then_initialize_filtered_thickness(self):
        with mock.patch.multiple(ReservoirGrid, from_csv=mock.DEFAULT, match_nodes=mock.DEFAULT) as mock_rgrid:
            thickness_grid = mock.MagicMock()
            mock_rgrid['from_csv'].return_value = thickness_grid
            ThinSheetModel(paths=dict(thickness='a', compressibility='b'), sigmas=dict(thickness=30))
        self.assertEqual(mock_rgrid['from_csv'].call_count, 2)
        thickness_grid.gaussian_filter.assert_called_once_with(30)
        mock_rgrid['match_nodes'].assert_called_once()

    def test_when_initializing_with_pressure_and_gaussian_filter_then_pressure_filter_is_adapted(self):
        pressure_grid = mock.MagicMock()
        pressure_grid.attributes = np.array([[[1, 1], [1, 1]], [[0, 0], [0, 0]]])
        with mock.patch.object(ReservoirGrid, 'match_nodes') as mock_rgrid:
            with mock.patch.object(ReservoirGridEpochs, 'from_csv') as mock_rgride:
                mock_rgride.return_value = pressure_grid
                with mock.patch.object(ReservoirGridEpochs, '__init__') as mock_rgride_init:
                    mock_rgride_init.return_value = None
                    ThinSheetModel(paths=dict(pressure='a'), sigmas=dict(pressure=30))
        mock_rgrid.assert_called_once()
        pressure_grid.snap.assert_called_once()
        pressure_grid.adapted_gaussian_filter.assert_called_once_with(30)
        mock_rgride.assert_called_once()
        mock_rgride_init.assert_called_once()

    def test_when_initializing_with_pressure_then_depletion_is_computed_from_first_epoch(self):
        pressure_grid = mock.MagicMock()
        pressure_grid.attributes = np.array([[[1, 1], [1, 1]], [[0, 0], [0, 0]]])
        with mock.patch.object(ReservoirGrid, 'match_nodes') as mock_rgrid:
            with mock.patch.object(ReservoirGridEpochs, 'from_csv') as mock_rgride:
                mock_rgride.return_value = pressure_grid
                with mock.patch.object(ReservoirGridEpochs, '__init__') as mock_rgride_init:
                    mock_rgride_init.return_value = None
                    ThinSheetModel(paths=dict(pressure='a'))
        mock_rgrid.assert_called_once()
        mock_rgride.assert_called_once()
        mock_rgride_init.assert_called_once()

        (_, _, _, act_attrs), _ = mock_rgride_init.call_args
        np.testing.assert_equal(act_attrs, np.array([[[0, 0], [0, 0]], [[.1, .1], [.1, .1]]]))

    def test_when_initializing_with_pressure_then_negative_depletion_is_set_to_zero(self):
        pressure_grid = mock.MagicMock()
        pressure_grid.attributes = np.array([[[1, 1], [1, 1]], [[0, 2], [0, 0]]])
        with mock.patch.object(ReservoirGrid, 'match_nodes') as mock_rgrid:
            with mock.patch.object(ReservoirGridEpochs, 'from_csv') as mock_rgride:
                mock_rgride.return_value = pressure_grid
                with mock.patch.object(ReservoirGridEpochs, '__init__') as mock_rgride_init:
                    mock_rgride_init.return_value = None
                    ThinSheetModel(paths=dict(pressure='a'))
        mock_rgrid.assert_called_once()
        mock_rgride.assert_called_once()
        mock_rgride_init.assert_called_once()

        (_, _, _, act_attrs), _ = mock_rgride_init.call_args
        np.testing.assert_equal(act_attrs, np.array([[[0, 0], [0, 0]], [[.1, 0], [.1, .1]]]))

    def test_when_initializing_with_pressure_then_depletion_has_same_nodes_and_names(self):
        xnodes, ynodes, names = np.array([0, 1]), np.array([1, 1]), ('a', 'b')
        pressure_grid = mock.MagicMock()
        pressure_grid.xnodes, pressure_grid.ynodes, pressure_grid.names = xnodes, ynodes, names
        pressure_grid.attributes = np.array([[[1, 1], [1, 1]], [[0, 0], [0, 0]]])
        with mock.patch.object(ReservoirGrid, 'match_nodes') as mock_rgrid:
            with mock.patch.object(ReservoirGridEpochs, 'from_csv') as mock_rgride:
                mock_rgride.return_value = pressure_grid
                with mock.patch.object(ReservoirGridEpochs, '__init__') as mock_rgride_init:
                    mock_rgride_init.return_value = None
                    ThinSheetModel(paths=dict(pressure='a'))
        mock_rgrid.assert_called_once()
        mock_rgride.assert_called_once()
        mock_rgride_init.assert_called_once()

        (act_xnodes, act_ynodes, act_names, _), _ = mock_rgride_init.call_args
        np.testing.assert_equal(act_xnodes, xnodes)
        np.testing.assert_equal(act_ynodes, ynodes)
        self.assertTupleEqual(act_names, names)


if __name__ == '__main__':
    unittest.main()
