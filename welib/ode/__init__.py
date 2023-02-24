# 


def integrate(f, t, x0, method='rk_nystrom', **kwargs):

    if method=='rk_nystrom':

        from .rk_nystrom import rk_nystrom_integrate
        return rk_nystrom_integrate(f, t, x0, acc0=None)

    elif method=='rk_dormand_prince':

        from .rk_dormand_prince import rk_dormand_prince_integrate
        return rk_dormand_prince_integrate(f, t, x0)

    else:
        raise NotImplementedError('Method {}'.format(method))

