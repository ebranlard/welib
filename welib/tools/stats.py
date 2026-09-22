""" 
Set of tools for statistics 
  - measures (R^2, RMSE)
  - pdf distributions
  - Binning 

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    from numpy import trapezoid
except:
    from numpy import trapz as trapezoid
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# --------------------------------------------------------------------------------}
# --- Stats measures 
# --------------------------------------------------------------------------------{
STAT_PRINT_CONFIG = { 
                 # Pretty name      # Latex fmt                                             # Plan format
    'sigratio' : (r'sigRatio'                    , r'$\sigma_\mathrm{{est}}/\sigma_\mathrm{{ref}} = {:.3f}$' , 'std ratio (est/ref)={:.3f}') , 
    'eps'      : (r'Rel. Err. $\epsilon$'        , r'$\epsilon={:.1f}\%$'                                , 'eps={:.1f}%')                , 
    'r2'       : (r'$R^2$'                       , r'$R^2={:.3f}$'                                       , 'R^2={:.3f}')                 , 
    'epsleq'   : (r'epsLeq'                      , r'$\epsilon L_{{eq}}={:.1f}\%$'                         , 'eps L_{eq}={:.1f}%')         , 
    'pearsonr' : (r'Pearson $\rho_{xy}(0)$'      , r'$\rho_{{xy}}(0)={:.3f}$'                                   , 'rho(0)={:.3f}')                   , 
    'spearmanr': (r'Spearman $\rho_S$'           , r'$\rho_S={:.3f}$'                                    , 'rho_S={:.3f}')                 , 
    'xcorr_max': (r'Max. Xcorr, $\rho_{xy,max}$' , r'$\rho_m={:.3f}$'                                    , 'xcorr={:.3f}')               , 
}

def comparison_stats(t1, y1, t2, y2, stats='sigRatio,eps,R2', method='mean', absVal=True, latex=True):
    """
    y1: ref
    y2: other

    """
    from welib.tools.fatigue import equivalent_load

    sp=stats.split(',')
    statsD = {}
    sStats=[]

    t1=np.asarray(t1).astype(float)
    y1=np.asarray(y1).astype(float)
    t2=np.asarray(t2).astype(float)
    y2=np.asarray(y2).astype(float)

    # Loop on statistics requested
    for stat_lab in sp:
        s= stat_lab.strip().lower()
        if s=='sigratio':
            # Ratio of standard deviation:
            sig_ref = float(np.nanstd(y1))
            sig_est = float(np.nanstd(y2))
            try:
                r_sig = sig_est/sig_ref
            except:
                r_sig = np.nan
            statsD[stat_lab] = r_sig

        elif s=='eps':
            # Mean relative error
            eps     = float(mean_rel_err(t1, y1, t2, y2, method=method, absVal=absVal))
            statsD[stat_lab] = eps

        elif s=='r2':
            # Rsquare
            R2 = float(rsquare(y_ref=y1, y_sim=y2)[0])
            statsD[stat_lab] = R2

        elif s=='epsleq':
            Leq1 = equivalent_load(t1, y1, m=5, bins=100, method='fatpack')
            Leq2 = equivalent_load(t2, y2, m=5, bins=100, method='fatpack')
            epsLeq = (Leq2-Leq1)/Leq1*100
            statsD[stat_lab] = epsLeq

        elif s in ['pearsonr', 'rho_xy(0)']:
            # Pearson is nothing more than the corss-correlation coefficient at zero lag
            try:
                r_val = float(pearsonr(y1, y2)[0])
            except Exception:
                r_val = np.nan
            statsD[stat_lab] = r_val # "r"

        elif s in ['spearmanr']:
            try:
                rho_val = float(spearmanr(y1, y2)[0])
            except Exception:
                rho_val = np.nan
            statsD[stat_lab] = rho_val # "rho"

        elif s in ['xcorr_max']:
            y1_norm = y1 - np.nanmean(y1)
            y2_norm = y2 - np.nanmean(y2)
            denom = np.nanstd(y1) * np.nanstd(y2) * len(y1)
            if denom != 0 and not np.isnan(denom):
                corr = np.correlate(y1_norm, y2_norm, mode='full') / denom
                xcorr_max = float(np.nanmax(corr))
            else:
                xcorr_max = np.nan
            statsD[stat_lab] = xcorr_max

        else:
            raise NotImplementedError()
        # --- 
        val = statsD[stat_lab]
        ss, latex_fmt, plain_fmt = STAT_PRINT_CONFIG[s]
        fmt = latex_fmt if latex else plain_fmt
        sStats += [fmt.format(val)]

    sStats=' - '.join(sStats)
    return statsD, sStats

def comparison_stats_grouped(sample_ID, y_ref, y_sim, mean_var=None, stats='sigRatio,eps,R2', method='mean', absVal=True, latex=True):
    """ Split samples based on sample_ID then computes average stats """
    sample_ID = np.asarray(sample_ID)
    y_ref = np.asarray(y_ref).astype(float)
    y_sim = np.asarray(y_sim).astype(float)

    unique_ids = np.unique(sample_ID)
    all_stats_dicts = []
    if mean_var is None:
        sample_ID

    # --- Split signal based on sample ID, and compute stats 
    mean_val =[]
    for uid in unique_ids:
        mask = (sample_ID == uid)
        sub_ref = y_ref[mask]
        sub_sim = y_sim[mask]
        sub_mean = mean_var[mask]
        sub_t = np.arange(len(sub_ref))
        mean_val.append(np.mean(sub_mean))

        s_dict, sStats = comparison_stats(sub_t, sub_ref, sub_t, sub_sim, stats=stats, method=method, absVal=absVal, latex=latex)
#         print('sStats', sStats)
        all_stats_dicts.append(s_dict)

    # --- Assemble values into a nice dict with numpy arrays
    dict_sample = {'ID': list(unique_ids)}
    keys = all_stats_dicts[0].keys() if all_stats_dicts else []
    for k in keys:
        dict_sample[k] = [d.get(k, np.nan) for d in all_stats_dicts]
    dict_sample['mean_val'] = mean_val

    # --- Take the average 
    mean_stats = {}
    for k in keys:
        vals = [d[k] for d in all_stats_dicts if k in d and not np.isnan(d[k])]
        mean_stats[k] = float(np.mean(vals)) if vals else np.nan

    # --- Write corredponding string
    sp = stats.split(',')
    sStats = []
    for s in sp:
        s = s.strip().lower()
        if s in STAT_PRINT_CONFIG:
            sss, latex_fmt, plain_fmt = STAT_PRINT_CONFIG[s]
            val = mean_stats.get(s, np.nan)
            fmt = latex_fmt if latex else plain_fmt
            sStats.append(fmt.format(val))
        else:
            raise NotImplementedError(s)

    sStats = ' - '.join(sStats)
    return mean_stats, sStats, dict_sample



def allclose_errors(actual, desired):
    actual = np.asarray(actual)
    desired = np.asarray(desired)

    # Absolute error element-wise
    abs_err = np.abs(actual - desired)
    max_abs_err = np.max(abs_err)

    # Relative error element-wise (handling division by zero safely)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_err = abs_err / np.abs(desired)
        # Filter out NaNs/Infs that occur where desired == 0
        rel_err_clean = np.where(np.isfinite(rel_err), rel_err, 0.0)
        max_rel_err = np.max(rel_err_clean)

    return max_abs_err, max_rel_err


def rsquare(y_ref, y_sim, c = True): 
    """ Compute coefficient of determination of data fit model and RMSE
    [r2 rmse] = rsquare(y_ref,y_sim)
    [r2 rmse] = rsquare(y_ref,y_sim,c)
    RSQUARE computes the coefficient of determination (R-square) value from
    actual data Y_REF and model data Y_SIM. The code uses a general version of
    R-square, based on comparing the variability of the estimation errors
    with the variability of the original values. RSQUARE also outputs the
    root mean squared error (RMSE) for the user's convenience.
    Note: RSQUARE ignores comparisons involving NaN values.
    INPUTS
      Y_REF     : Actual data
      Y_SIM     : Model fit
    
    # OPTION
      C         : Constant term in model
                    R-square may be a questionable measure of fit when no
                  constant term is included in the model.
      [DEFAULT] TRUE : Use traditional R-square computation
                FALSE : Uses alternate R-square computation for model
                       without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
    # OUTPUT
      R2        : Coefficient of determination
      RMSE      : Root mean squared error """
    # Sanity
    if not np.all(y_ref.shape == y_sim.shape) :
        raise Exception('Y_REF and Y_SIM must be the same size')
    y_ref = np.asarray(y_ref).astype(float)
    y_sim = np.asarray(y_sim).astype(float)
    # Check for NaN
    tmp = np.logical_not(np.logical_or(np.isnan(y_ref),np.isnan(y_sim))) 
    y_ref = y_ref[tmp]
    y_sim = y_sim[tmp]
    if c:
        denom = np.sum((y_ref-np.mean(y_ref))** 2)
        if abs(denom)>0:
            r2 = max(0,1-np.sum((y_ref-y_sim)**2)/denom)
        else:
            r2 = np.inf
    else:
        denom = np.sum((y_ref) ** 2)
        if abs(denom)>0:
            r2 = 1 - np.sum((y_ref - y_sim) ** 2) /denom
        else:
            r2 = np.inf
        if r2 < 0:
            import warnings
            warnings.warn('Consider adding a constant term to your model')
            r2 = 0
    rmse = np.sqrt(np.mean((y_ref - y_sim) ** 2))
    return r2,rmse





def mean_rel_err(t1=None, y1=None, t2=None, y2=None, method='meanabs', verbose=False, varname='', absVal=True, tRange=None):
    """ 
    return mean relative error in % 

    Methods: 
      'mean'   : 100 * |y1-y2|/mean(y1)
      'meanabs': 100 * |y1-y2|/mean(|y1|)
      'minmax': y1 and y2 scaled between 0.5 and 1.5
                |y1s-y2s|/|y1|
      '0-2': signals are scalled between 0 & 2
    """
    def myabs(y):
        if absVal:
            return np.abs(y)
        else:
            return y

    if tRange is not None and t1 is not None:
        b = np.logical_and(t1>tRange[0], t1<tRange[1])
        if sum(b)>0:
            t1 = t1[b]
            y1 = y1[b]
            if t2 is not None:
                b=np.logical_and(t2>tRange[0], t2<tRange[1])
                t2 = t2[b]
                y2 = y2[b]


    if t1 is None and t2 is None:
        pass
    else:
        if len(y1)!=len(y2):
            y2=np.interp(t1, t2, y2)


#     print('Mean rel error {:7.2f} %'.format( meanrelerr))
#     return meanrelerr,meanrelerr0

    if method=='mean':
        # Method 1 relative to mean
        ref_val = np.nanmean(y1)
        if abs(ref_val)>0:
            meanrelerr = np.nanmean(myabs(y2-y1)/ref_val)*100 
        else:
            meanrelerr = np.nan
    elif method=='meanabs':
        ref_val = np.nanmean(abs(y1))
        if abs(ref_val)>0:
            meanrelerr = np.nanmean(myabs(y2-y1)/ref_val)*100 
        else:
            meanrelerr = np.nan
    elif method=='loc':
        meanrelerr = np.nanmean(myabs(y2-y1)/abs(y1))*100 
    elif method=='minmax':
        # Method 2 scaling signals
        Min=min(np.nanmin(y1), np.nanmin(y2))
        Max=max(np.nanmax(y1), np.nanmax(y2))
        y1=(y1-Min)/(Max-Min)+0.5
        y2=(y2-Min)/(Max-Min)+0.5
        meanrelerr = np.nanmean(myabs(y2-y1)/np.abs(y1))*100 
    elif method=='1-2':
        # transform values from 1 to 2
        Min=min(np.nanmin(y1), np.nanmin(y2))
        Max=max(np.nanmax(y1), np.nanmax(y2))
        if Max==Min:
            Max=Min+1
        y1 = (y1-Min)/(Max-Min)+1
        y2 = (y2-Min)/(Max-Min)+1
        meanrelerr = np.nanmean(myabs(y2-y1)/np.abs(y1))*100

    else:
        raise Exception('Unknown method',method)

    if verbose:
        if len(varname)>0:
            print('Mean rel error {:15s} {:7.2f} %'.format(varname, meanrelerr))
        else:
            print('Mean rel error {:7.2f} %'.format( meanrelerr))
    return meanrelerr


# --------------------------------------------------------------------------------}
# --- PDF 
# --------------------------------------------------------------------------------{
def pdf(y, method='histogram', n=50, **kwargs):
    """ 
    Compute the probability density function.
    Wrapper over the different methods present in this package
    """
    if method =='sns':
        xh, yh = pdf_sns(y, nBins=n, **kwargs)
    elif method =='gaussian_kde':
        xh, yh = pdf_gaussian_kde(y, nOut=n, **kwargs)
    elif method =='histogram':
        xh, yh = pdf_histogram(y, nBins=n, **kwargs)
    else:
        raise NotImplementedError(f'pdf method: {method}')
    return xh, yh


def pdf_histogram(y,nBins=50, norm=True, count=False):
    yh, xh = np.histogram(y[~np.isnan(y)], bins=nBins)
    dx   = xh[1] - xh[0]
    xh  = xh[:-1] + dx/2
    if count:
        yh  = yh / (len(n)*dx) # TODO DEBUG /VERIFY THIS
    else:
        yh  = yh / (nBins*dx) 
    if norm:
        try:
            yh=yh/trapezoid(yh,xh)
        except:
            yh=yh/np.trapz(yh,xh)
    return xh,yh

def pdf_gaussian_kde(data, bw='scott', nOut=100, cut=3, clip=(-np.inf,np.inf)):
    """ 
    Returns a smooth probability density function (univariate kernel density estimate - kde) 
    Inspired from `_univariate_kdeplot` from `seaborn.distributions`

    INPUTS:
        bw:  float defining bandwidth or method (string) to find it (more or less sigma)   
        cut: number of bandwidth kept for x axis (e.g. 3 sigmas)
        clip: (xmin, xmax) values
    OUTPUTS:
        x, y: where y(x) = pdf(data)
    """
    from scipy import stats
    from six import string_types

    data = np.asarray(data)
    data = data[~np.isnan(data)]
    # Gaussian kde
    kde  = stats.gaussian_kde(data, bw_method = bw)
    # Finding a relevant support (i.e. x values)
    if isinstance(bw, string_types):
        bw_ = "scotts" if bw == "scott" else bw
        bw = getattr(kde, "%s_factor" % bw_)() * np.std(data)
    x_min = max(data.min() - bw * cut, clip[0])
    x_max = min(data.max() + bw * cut, clip[1])
    x = np.linspace(x_min, x_max, nOut)
    # Computing kde on support
    y = kde(x)
    return x, y


def pdf_sklearn(y):
    #from sklearn.neighbors import KernelDensity
    #kde = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(y) #you can supply a bandwidth
    #x=np.linspace(0,5,100)[:, np.newaxis]
    #log_density_values=kde.score_samples(x)
    #density=np.exp(log_density)
    pass

def pdf_sns(y,nBins=50):
    import seaborn.apionly as sns
    hh=sns.distplot(y,hist=True,norm_hist=False).get_lines()[0].get_data()
    xh=hh[0]
    yh=hh[1]
    return xh,yh

# --------------------------------------------------------------------------------}
# --- Binning 
# --------------------------------------------------------------------------------{
def bin_DF(df, xbins, colBin, stats=None):
    """ 
    Perform bin averaging of a dataframe
    INPUTS:
      - df   : pandas dataframe
      - xBins: end points delimiting the bins, array of ascending x values
      - colBin: column name (string) of the dataframe, used for binning 
    OUTPUTS:
       binned dataframe, with additional columns 'Counts' for the number 

    """
    if stats is None:
        stats=['avg']
    if not isinstance(stats, list):
        stats=[stats]
    if colBin not in df.columns.values:
        raise Exception('The column `{}` does not appear to be in the dataframe'.format(colBin))
    xmid      = (xbins[:-1]+xbins[1:])/2
    df['Bin'] = pd.cut(df[colBin], bins=xbins, labels=xmid ) # Adding a column that has bin attribute
    dfs=[]
    df3  = df.groupby('Bin', observed=False)
    for stat in stats:
        if stat=='avg' or stat=='mean':
            df2  = df3.mean()  # mean by bin
        elif stat=='std':
            df2  = df3.std()   # std by bin
        elif stat=='min':
            df2  = df3.min()   # min by bin
        elif stat=='max':
            df2  = df3.max()   # min by bin
        else:
            raise NotImplementedError(f'Stat {stat}')
        df2  = df2.reindex(xmid) # Just in case some bins are missing (will be nan)
        dfs.append(df2)
    # Adding counts to first df
    df['Counts'] = 1
    dfCount=df[['Counts','Bin']].groupby('Bin', observed=False).sum()
    dfs[0]['Counts'] = dfCount['Counts']
    return dfs



def bin_signal(x, y, xbins=None, stats=None, nBins=None):
    """ 
    Perform bin averaging of a signal
    INPUTS:
      - x: x-values 
      - y: y-values, signal values
      - xBins: end points delimiting the bins, array of ascending x values
    OUTPUTS:
      - xBinned, yBinned

    """
    if stats is None:
        stats=['avg']
    if not isinstance(stats, list):
        stats=[stats]
    if xbins is None:
        xmin, xmax = np.min(x), np.max(x)
        dx = (xmax-xmin)/nBins
        xbins=np.arange(xmin, xmax+dx/2, dx)
    df = pd.DataFrame(data=np.column_stack((x,y)), columns=['x','y'])
    dfs = bin_DF(df, xbins, colBin='x', stats=stats)
    if len(stats)>1:
        raise NotImplementedError('bin_signal for multiple stats')
    else:
        return dfs[0]['x'].values, dfs[0]['y'].values



def bin2d_signal(x, y, z, xbins=None, ybins=None, nXBins=None, nYBins=None):
    """ 
    Bin signal z based on x and y values using xbins and ybins

    """
    if xbins is None:
        xmin, xmax = np.min(x), np.max(x)
        dx = (xmax-xmin)/nXBins
        xbins=np.arange(xmin, xmax+dx/2, dx)
    if ybins is None:
        ymin, ymax = np.min(y), np.max(y)
        dy = (ymax-ymin)/nYBins
        ybins=np.arange(ymin, ymax+dy/2, dy)

    x = np.asarray(x).flatten()
    y = np.asarray(y).flatten()
    z = np.asarray(z).flatten()

    Counts = np.zeros((len(xbins)-1, len(ybins)-1))
    XMean  = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan
    YMean  = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan
    ZMean  = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan
    ZStd   = np.zeros((len(xbins)-1, len(ybins)-1))*np.nan

    xmid = xbins[:-1] + np.diff(xbins)/2
    ymid = ybins[:-1] + np.diff(ybins)/2
    YMid, XMid = np.meshgrid(ymid, xmid)

    for ixb, xb in enumerate(xbins[:-1]):
        print(ixb)
        bX = np.logical_and(x >= xb, x <= xbins[ixb+1]) # TODO decide on bounds
        for iyb, yb in enumerate(ybins[:-1]):
            bY = np.logical_and(y >= yb, y <= ybins[iyb+1]) # TODO decide on bounds

            bXY = np.logical_and(bX, bY)
            Counts[ixb, iyb] = sum(bXY)
            if Counts[ixb,iyb]>0:
                ZMean [ixb, iyb] = np.mean(z[bXY])
                ZStd  [ixb, iyb] = np.std( z[bXY])
                XMean [ixb, iyb] = np.mean(x[bXY])
                YMean [ixb, iyb] = np.mean(y[bXY])

    return XMean, YMean, ZMean, ZStd, Counts, XMid, YMid






def azimuthal_average_DF(df, psiBin=np.arange(0,360+1,10), colPsi='Azimuth_[deg]', tStart=None, colTime='Time_[s]'):
    """ 
    Average a dataframe based on azimuthal value
    Returns a dataframe with same amount of columns as input, and azimuthal values as index
    """
    if tStart is not None:
        if colTime not in df.columns.values:
            raise Exception('The column `{}` does not appear to be in the dataframe'.format(colTime))
        df=df[ df[colTime]>tStart].copy()

    dfPsi= bin_DF(df, psiBin, colPsi, stats=['avg'])[0]
    if np.any(dfPsi['Counts']<1):
        print('[WARN] some bins have no data! Increase the bin size.')

    return dfPsi


def azimuthal_std_DF(df, psiBin=np.arange(0,360+1,10), colPsi='Azimuth_[deg]', tStart=None, colTime='Time_[s]'):
    """ 
    Average a dataframe based on azimuthal value
    Returns a dataframe with same amount of columns as input, and azimuthal values as index
    """
    if tStart is not None:
        if colTime not in df.columns.values:
            raise Exception('The column `{}` does not appear to be in the dataframe'.format(colTime))
        df=df[ df[colTime]>tStart].copy()

    dfPsi= bin_DF(df, psiBin, colPsi, stats=['std'])[0]
    if np.any(dfPsi['Counts']<1):
        print('[WARN] some bins have no data! Increase the bin size.')

    return dfPsi




def plot_yy(y_sim, y_ref, ax=None, 
            label = None, sc_label=None, bin_label=None, 
            scatter=True, sc_color=None, sc_size=None,  # Scatter Options
            sc_alpha=0.5, sc_marker='o'                 ,  # Scatter Options
            nBins=20, bin_color=None, bin_ls='-', bin_marker='o', # Bin Options
            bin_markeredgecolor='k', bin_markersize=None,                        # Bin Options
            pdf_x = False,                                             # PDF options
            lims=None,
            lg=None, lg_fs=11, lg_statsInBox=True, lg_loc='left',
            stats='eps,R2'):
    """ Perform a y-y plot, with"""
    
    if ax is None:
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        
    y_sim = np.asarray(y_sim)
    y_ref = np.asarray(y_ref)


    # --- y-y black reference
    all_vals = np.concatenate([y_ref, y_sim])
    valid_vals = all_vals[~np.isnan(all_vals)]
    if len(valid_vals) > 0:
        if lims is None:
            lims = [np.min(valid_vals)*1.05, np.max(valid_vals)*1.05]
            ax.set_xlim(lims)
            ax.set_ylim(lims)
        ax.plot(lims, lims, color='black', linestyle='--', alpha=0.7)
    
    # --- Compute stats 
    sStats = ""
    if stats:
        t = np.arange(len(y_ref))
        _, sStats = comparison_stats(t, y_ref, t, y_sim, stats=stats)

    # --- Determine where to put stats (in labels or box)
    if stats and (not lg_statsInBox):
        if sc_label is None and bin_label is None:
            sc_label = sStats if scatter else None
            bin_label = sStats if (nBins and scatter is None) else None
        elif sc_label is not None:
            sc_label += ' - '+sStats 
        elif bin_label is not None:
            bin_label += ' - '+sStats 

    # --- Scatter plot
    if scatter:
        ax.scatter(y_ref, y_sim, c=sc_color, s=sc_size, label=sc_label, alpha=sc_alpha, marker=sc_marker)
        
    # --- Bin plot
    if nBins is not None:
        xBinned, yBinned = bin_signal(y_ref, y_sim, nBins=nBins, stats=['avg'])
        ax.plot(xBinned, yBinned, color=bin_color, linestyle=bin_ls, marker=bin_marker,
                markeredgecolor=bin_markeredgecolor, markersize=bin_markersize,
                label=bin_label)
        
        
    # -- Esthetics
    ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
    if lg:
        if bin_label or sc_label:
            ax.legend(loc=lg_loc, fontsize=lg_fs)
    if lg_statsInBox:
        # Split the combined stats string and join each item with a newline character
        if label is not None:
            sStats = label + ' - ' + sStats

        lines = sStats.split(' - ')
        if lines:
            if lines[0].startswith('$') and lines[0].endswith('$'):
                inner = lines[0][1:-1]
            else:
                inner = lines[0]
            lines[0] = f"$\\mathbf{{{inner}}}$"
        stats_text = '\n'.join(lines)

        # Place the text box inside your plot axes
        ax.text(
            0.05, 0.95, stats_text,
            transform=ax.transAxes,
            fontsize=lg_fs,
            verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray')
        )
        
    return ax

if __name__ == '__main__':
    # Dummy test execution
    np.random.seed(42)
    y_ref = np.linspace(0, 10, 100)
    y_sim = y_ref + np.random.normal(0, 1, 100)

    fig, ax = plt.subplots(figsize=(6, 6))
    plot_yy(y_sim, y_ref, ax=ax, scatter=True, nBins=10)
    plt.xlabel("Reference Values")
    plt.ylabel("Simulation Values")
    plt.title("Dummy Test Example")
    plt.show()
