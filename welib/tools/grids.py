import numpy as np
import matplotlib.pyplot as plt

class BoundingBox:
    def __init__(self, x, y, z=None):
        """Initialize from arrays of coordinates x,y,(z)."""
        self.is3D = z is not None
        self.xmin, self.xmax = np.min(x), np.max(x)
        self.ymin, self.ymax = np.min(y), np.max(y)
        if self.is3D:
            self.zmin, self.zmax = np.min(z), np.max(z)

    def lines(self):
        """Return line segments (arrays of shape (N,2) or (N,3))
        for plotting the bounding box edges."""
        if not self.is3D:
            # 2D rectangle (closed loop)
            pts = np.array([
                [self.xmin, self.ymin],
                [self.xmax, self.ymin],
                [self.xmax, self.ymax],
                [self.xmin, self.ymax],
                [self.xmin, self.ymin]
            ])
            return pts
        else:
            # 3D: 8 corners
            corners = np.array([
                [self.xmin, self.ymin, self.zmin],
                [self.xmax, self.ymin, self.zmin],
                [self.xmax, self.ymax, self.zmin],
                [self.xmin, self.ymax, self.zmin],
                [self.xmin, self.ymin, self.zmax],
                [self.xmax, self.ymin, self.zmax],
                [self.xmax, self.ymax, self.zmax],
                [self.xmin, self.ymax, self.zmax]
            ])
            # faces as loops of 4 corners (closed)
            faces = [
                [0,1,2,3,0], # bottom
                [4,5,6,7,4], # top
                [0,1,5,4,0], # front
                [1,2,6,5,1], # right
                [2,3,7,6,2], # back
                [3,0,4,7,3]  # left
            ]

            X, Y, Z = [], [], []
            for f in faces:
                X.extend(corners[f,0]); X.append(np.nan)  # NaN breaks the line
                Y.extend(corners[f,1]); Y.append(np.nan)
                Z.extend(corners[f,2]); Z.append(np.nan)

            points = np.zeros((len(X), 3) )
            points[:,0] = X
            points[:,1] = Y
            points[:,2] = Z
            return points

    def contains(self, other, strict: bool = False):
        """Return True if self contains other bounding box.

        Parameters
        ----------
        other : BoundingBox
            The other bounding box to check.
        strict : bool, optional
            If True, requires strict containment (no shared boundaries).
            Default is False.
        """
        if self.is3D != other.is3D:
            raise Exception('Cannot compare a 2D and 3D box')
        if strict:
            ok = (self.xmin < other.xmin) and (self.xmax > other.xmax) \
                 and (self.ymin < other.ymin) and (self.ymax > other.ymax)
            if self.is3D:
                ok = ok and (self.zmin < other.zmin) and (self.zmax > other.zmax)
        else:
            ok = (self.xmin <= other.xmin) and (self.xmax >= other.xmax) \
                 and (self.ymin <= other.ymin) and (self.ymax >= other.ymax)
            if self.is3D:
                ok = ok and (self.zmin <= other.zmin) and (self.zmax >= other.zmax)
        return ok

    def plot(self, ax, plane=None, **kwargs):
        """Plot bounding box in the chosen plane."""
        pts = self.lines()
        if plane is None:
            if not self.is3D:
                ax.plot(pts[:,0], pts[:,1], **kwargs)
            else:
                ax.plot(pts[:,0], pts[:,1], pts[:,2], **kwargs)
            return
        if plane == "XY":
            X, Y = pts[:,0],pts[:,1]
        elif plane == "XZ" and self.is3D:
            X, Y = pts[:,0],pts[:,2]
        elif plane == "YZ" and self.is3D:
            X, Y = pts[:,1],pts[:,2]
        else:
            raise ValueError(f"Invalid plane {plane} for this bounding box")
        ax.plot(X, Y, **kwargs)



class RegularGrid:
    def __init__(self, x0, nx, dx, y0=None, ny=None, dy=None, z0=None, nz=None, dz=None):
        """
        Define a regular grid in 2D or 3D.

        Parameters
        ----------
        x0, nx, dx : float, int, float
            Origin, number of points, spacing along x.
        y0, ny, dy : optional
            Same for y.
        z0, nz, dz : optional
            Same for z.
        """
        self.x0, self.nx, self.dx = x0, nx, dx
        self.y0, self.ny, self.dy = y0, ny, dy
        self.z0, self.nz, self.dz = z0, nz, dz
        self.is3D = z0 is not None and nz is not None and dz is not None

        # Coordinates
        self.x = x0 + np.arange(nx) * dx
        if y0 is not None and ny is not None and dy is not None:
            self.y = y0 + np.arange(ny) * dy
        else:
            self.y = None
        if self.is3D:
            self.z = z0 + np.arange(nz) * dz
        else:
            self.z = None

        # Bounding box
        if self.is3D:
            self.bb = BoundingBox(self.x, self.y, self.z)
        else:
            self.bb = BoundingBox(self.x, self.y)

    def contains_grid(self, other, strict=False):
        """Check if this grid fully contains another grid."""
        return self.bb.contains(other.bb, strict=strict)

    def contains_bb(self, other_bb, strict=False):
        """Check if this grid contains a given bounding box."""
        return self.bb.contains(other_bb, strict=strict)

    def contains_p(self, p):
        """Check if this grid contains a point p=(x,y) or (x,y,z)."""
        if self.is3D and len(p) == 3:
            x, y, z = p
            return (self.bb.xmin <= x <= self.bb.xmax and
                    self.bb.ymin <= y <= self.bb.ymax and
                    self.bb.zmin <= z <= self.bb.zmax)
        elif not self.is3D and len(p) == 2:
            x, y = p
            return (self.bb.xmin <= x <= self.bb.xmax and
                    self.bb.ymin <= y <= self.bb.ymax)
        else:
            return False

    def plot(self, ax, plane="XY", grid=True, color=(0.3,0.3,0.3), optsGd=None, optsBB=None):
        """Plot grid lines and bounding box projection in a given plane."""
        optsGdLoc=dict(ls='-', color=color, lw=0.3)
        optsGd=optsGdLoc if optsGd is None else optsGdLoc.update(optsGd)
        optsBBLoc=dict(ls='-', color=color, lw=1.0)
        optsBB=optsBBLoc if optsBB is None else optsBBLoc.update(optsBB)

        if plane == "XY":
            X, Y = self.x, self.y
        elif plane == "XZ" and self.is3D:
            X, Y = self.x, self.z
        elif plane == "YZ" and self.is3D:
            X, Y = self.y, self.z
        else:
            raise ValueError(f"Invalid plane {plane} for this grid")

        # Grid lines
        if grid:
            ax.vlines(X, ymin=Y[0], ymax=Y[-1], **optsGdLoc)
            ax.hlines(Y, xmin=X[0], xmax=X[-1], **optsGdLoc)
        # Bounding box
        self.bb.plot(ax, plane=plane, **optsBBLoc)


if __name__ == "__main__":
    # --- 2D test ---
    x = np.random.randn(20)
    y = np.random.randn(20)
#     bb2d = BoundingBox(x, y)
# 
#     fig, ax = plt.subplots()
#     ax.scatter(x, y, label="points")
#     bb2d.plot(ax, color='r', lw=2, label="bounding box")
#     ax.legend()
#     plt.title("2D BoundingBox test")
#     plt.show()

    # --- 3D test ---
    z = np.random.randn(20)
    bb3d = BoundingBox(x, y, z)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='b', marker='o')
    bb3d.plot(ax, color='g', lw=2)
    plt.title("3D BoundingBox test")
    plt.show()

    # --- containment test ---
    x_small = np.random.uniform(-0.5, 0.5, 10)
    y_small = np.random.uniform(-0.5, 0.5, 10)
    z_small = np.random.uniform(-0.5, 0.5, 10)
    small_box = BoundingBox(x_small, y_small, z_small)
    print("bb3d contains small_box?", bb3d.contains(small_box))

