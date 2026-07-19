import matplotlib.pyplot as plt
from welib.airfoils.shapes import AirfoilShape
from welib.airfoils.naca import naca_shape

def write_coords(digits, n=600):
    x, y = naca_shape(digits, chord=1, n=n, sharp=True)
    arf = AirfoilShape(x=x, y=y, name='Naca'+digits)
    #print('x',arf.x)
    #print('y',arf.y)
    arf.write('_Naca{}.csv'.format(digits), format='csv', delim=' ')
    arf.plot_surfaces(title = 'Airfoils - Shapes - NACA '+ digits)
    #arf.plot()


if __name__ == '__main__':
    write_coords(digits='0012', n=600)
    write_coords(digits='0018', n=600)
    write_coords(digits='4403', n=600)
    write_coords(digits='4412', n=600)
    write_coords(digits='4416', n=600)
    write_coords(digits='4418', n=600)
    write_coords(digits='4430', n=600)


    plt.show()
if __name__ == '__test__':
    write_coords(digits='0018', n=36)
    pass
if __name__=="__export__":
    write_coords(digits='0018', n=36)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

