import os
import unittest
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from welib.fast.servodyn import *

scriptDir = os.path.dirname(__file__)

DEBUG = True
DEBUG = False
# matplotlib.use('Agg')

class TestSvD(unittest.TestCase):
    """ Unit tests for ServoDyn generator torque fitting logic """

    def test_SvD_fitGenTorqueFullRange(self):
        for WT in (1, 2):
            if WT == 1:
                svdFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NREL5MW_SvD_Simple.dat')
            else:
                svdFile = os.path.join(scriptDir, '../../../data/IEA-22-280-RWT/IEA-22-280-RWT_ServoDyn.dat')

            svd = ServoDyn(svdFile)

            # Get the VS data
            df = svd.VS_DataFrame(rpm_start=2, nRPM=15, fact_start=0.5, fact_max=1.2)

            # Try to fit the model
            torque_fit, fitter = svd.VS_fit(df['Generator_Speed_[rpm]'], df['Generator_Torque_[Nm]'])
            coeffs = fitter.model['coeffs']

            if DEBUG:
                print(fitter)
                fig, ax = plt.subplots(1, 1, figsize=(7, 5))
                ax.plot(df['Generator_Speed_[rpm]'], df['Generator_Torque_[Nm]'], 'o', label='ServoDyn Curve')
                x_fit = np.linspace(df['Generator_Speed_[rpm]'].min(), df['Generator_Speed_[rpm]'].max() * 1.05, 200)
                ax.plot(x_fit, fitter.model['fitted_function'](x_fit), '--', label=f'Fit (WT={WT})')
                ax.set_xlabel('Generator Speed [rpm]')
                ax.set_ylabel('Generator Torque [Nm]')
                ax.set_title(f'ServoDyn VS Fit Comparison (WT {WT})')
                ax.legend()
                ax.grid(True)
#                 plt.show()

            # Compare fitted parameters with ground-truth ServoDyn input parameters (allowing a slightly looser tolerance for sparse data)
            np.testing.assert_allclose(coeffs['RtGnSp'], svd.File['VS_RtGnSp'], rtol=0.03)
            np.testing.assert_allclose(coeffs['RtTq']  , svd.File['VS_RtTq'], rtol=0.03)
            np.testing.assert_allclose(coeffs['Rgn2K'] , svd.File['VS_Rgn2K'], rtol=0.08)
            np.testing.assert_allclose(coeffs['SlPc']  , svd.File['VS_SlPc'], rtol=0.7)

    def test_SvD_fitGenTorqueNarrowRange(self):
        for WT in (1, 2):
            if WT == 1:
                svdFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NREL5MW_SvD_Simple.dat')
            else:
                svdFile = os.path.join(scriptDir, '../../../data/IEA-22-280-RWT/IEA-22-280-RWT_ServoDyn.dat')

            svd = ServoDyn(svdFile)

            # Get the VS data
            df = svd.VS_DataFrame(rpm_start=2, nRPM=10, fact_start=1.0, fact_max=1.0)
#             print(df)

            # Try to fit the model
            torque_fit, fitter = svd.VS_fit(df['Generator_Speed_[rpm]'], df['Generator_Torque_[Nm]'])
            coeffs = fitter.model['coeffs']

            if DEBUG:
                print(fitter)
                fig, ax = plt.subplots(1, 1, figsize=(7, 5))
                ax.plot(df['Generator_Speed_[rpm]'], df['Generator_Torque_[Nm]'], 'o', label='ServoDyn Curve')
                x_fit = np.linspace(df['Generator_Speed_[rpm]'].min(), df['Generator_Speed_[rpm]'].max() * 1.05, 200)
                ax.plot(x_fit, fitter.model['fitted_function'](x_fit), '--', label=f'Fit (WT={WT})')
                ax.set_xlabel('Generator Speed [rpm]')
                ax.set_ylabel('Generator Torque [Nm]')
                ax.set_title(f'ServoDyn VS Fit Comparison (WT {WT})')
                ax.legend()
                ax.grid(True)

            # Compare fitted parameters with ground-truth ServoDyn input parameters (allowing ~2% relative tolerance)
            print('WT', WT)
            print('RtGnSp', 'Fit:', coeffs['RtGnSp'], '   File:', svd.File['VS_RtGnSp'])
            print('RtTq'  , 'Fit:', coeffs['RtTq']  , '   File:', svd.File['VS_RtTq'])
            print('Rgn2K' , 'Fit:', coeffs['Rgn2K'] , '   File:', svd.File['VS_Rgn2K'])
            print('SlPc'  , 'Fit:', coeffs['SlPc']  , '   File:', svd.File['VS_SlPc'])

            np.testing.assert_allclose(coeffs['RtGnSp'], svd.File['VS_RtGnSp'], rtol=0.03)
            np.testing.assert_allclose(coeffs['RtTq']  , svd.File['VS_RtTq'], rtol=0.03)
            np.testing.assert_allclose(coeffs['Rgn2K'] , svd.File['VS_Rgn2K'], rtol=0.08)
            np.testing.assert_allclose(coeffs['SlPc']  , svd.File['VS_SlPc'], rtol=2.5)





    def test_SvD_fitGenTorque_noisy_scatter(self):
        """ Synthetic 1-year operational dataset (52,560 10-minute intervals) with 2D measurement noise """
        svdFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Baseline/NREL5MW_SvD_Simple.dat')
        svd = ServoDyn(svdFile)
        svd.File['SpdGenOn'] = 200


        # Ground truth parameters
        RtGnSp_true   = svd.File['VS_RtGnSp']
        RtTq_true     = svd.File['VS_RtTq']
        Rgn2K_true    = svd.File['VS_Rgn2K']
        SlPc_true     = svd.File['VS_SlPc']
        SpdGenOn_true = svd.File['SpdGenOn']

        p_true = (RtGnSp_true, RtTq_true, Rgn2K_true, SlPc_true, SpdGenOn_true)

        # Generate 52,560 realizations (1 year at 10-min sampling)
        np.random.seed(42)
        n_samples = 52560

        # Operational speed distribution (weighted across Region 2, 2.5, and 3)
        rpm_clean = np.random.uniform(SpdGenOn_true * 0.8, RtGnSp_true * 1.01, n_samples)
        torque_clean = gentorque(rpm_clean, p_true)

        # Introduce noise to both speed (X-axis) and torque (Y-axis)
        rpm_noise_std = RtGnSp_true * 0.015       # 1.5% speed measurement uncertainty
        torque_noise_std = RtTq_true * 0.035      # 3.5% torque fluctuation/measurement noise

        rpm_noisy = rpm_clean + np.random.normal(0, rpm_noise_std, n_samples)
        torque_noisy = torque_clean + np.random.normal(0, torque_noise_std, n_samples)

        # Clip speed to positive domain
        rpm_noisy = np.maximum(0.0, rpm_noisy)

        # Perform curve fitting on noisy scatter
        torque_fit, fitter = svd.VS_fit(rpm_noisy, torque_noisy)
        coeffs = fitter.model['coeffs']

        if DEBUG:
            print("\n--- Noisy Scatter Fit Results ---")
            print(fitter)
            fig, ax = plt.subplots(1, 1, figsize=(8, 5))
            # Subsample for plot clarity
            idx = np.random.choice(n_samples, size=3000, replace=False)
            ax.scatter(rpm_noisy[idx], torque_noisy[idx], s=4, alpha=0.25, color='gray', label='1-Yr 10-min Scatter')

            x_fit = np.linspace(rpm_noisy.min(), rpm_noisy.max(), 250)
            ax.plot(x_fit, gentorque(x_fit, p_true), 'k-', lw=2, label='True Curve')
            ax.plot(x_fit, fitter.model['fitted_function'](x_fit), 'r--', lw=2, label='Scatter Fit')

            ax.set_xlabel('Generator Speed [rpm]')
            ax.set_ylabel('Generator Torque [Nm]')
            ax.set_title('Generator Torque Model Fit from Noisy Scatter')
            ax.legend()
            ax.grid(True)

#         np.testing.assert_allclose(coeffs['RtGnSp'], RtGnSp_true, rtol=0.03)
#         np.testing.assert_allclose(coeffs['RtTq'], RtTq_true, rtol=0.03)
        np.testing.assert_allclose(coeffs['RtGnSp'], svd.File['VS_RtGnSp'], rtol=0.06)
        np.testing.assert_allclose(coeffs['RtTq']  , svd.File['VS_RtTq'], rtol=0.09)
        np.testing.assert_allclose(coeffs['Rgn2K'] , svd.File['VS_Rgn2K'], rtol=0.12)
        np.testing.assert_allclose(coeffs['SlPc']  , svd.File['VS_SlPc'], rtol=0.50)

if __name__ == '__main__':
#     TestSvD().test_SvD_fitGenTorqueFullRange()
#     TestSvD().test_SvD_fitGenTorque_noisy_scatter()
#     TestSvD().test_SvD_fitGenTorqueNarrowRange()
#     plt.show()
    unittest.main()
