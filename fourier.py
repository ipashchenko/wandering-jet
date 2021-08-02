import numpy as np
import finufftpy


class FINUFFT:
    def __init__(self, uv, pixsize):
        self.pixsize = pixsize
        self.u = uv[:, 0]*pixsize*2*np.pi
        self.v = uv[:, 1]*pixsize*2*np.pi
        self.size_nu = len(uv)

    def forward(self, img):
        """
        :param img:
            2D array with shape (ms, mt)
        :return:
            c[j] =   SUM   img[k1,k2] exp(+/-i (k1 x[j] + k2 y[j])),  for j = 0,...,nj-1
                k1,k2

            where sum is over -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2
            and x and y are u & v.

        """
        vis = np.empty((self.size_nu,), dtype=np.complex128, order='F')
        img = np.array(img, dtype=np.complex128, order='F')
        code = finufftpy.nufft2d2(self.u, self.v, vis, isign=-1, eps=1e-8,
                                  f=img, fftw=1, upsampfac=1.25)
        return vis


class FINUFFT_NUNU:
    def __init__(self, uv, x, y):
        """

               nj-1
        f[k] = SUM c[j]*exp(+/-i*s[k]*x[j]+t[k]*y[j]),  for k = 0,...,nk-1
               j=0

        :param uv:
            Iterable of float uv-coordinates (length ``nk``)
        :param x:
            Iterable of float x-coordinates (length ``nj``). Multiplied on
            ``u`` in FT.
        :param y:
            Iterable of float y-coordinates (length ``nj``). Multiplied on
            ``v`` in FT.

        """

        u = uv[:, 0]*2*np.pi
        v = uv[:, 1]*2*np.pi
        UV_max = max(max(abs(u)), max(abs(v)))
        self.u = u/UV_max
        self.v = v/UV_max
        self.x = x*UV_max
        self.y = y*UV_max
        self.size_nu = len(uv)

    def forward(self, img):
        """
        :param img:
            Iterable of image pixel values (length ``nj``)
        :return:

                     nj-1
            f[k]  =  SUM   c[j] exp(+-i s[k] x[j] + t[k] y[j]),  for k = 0,...,nk-1
                     j=0

          Args:
            x     (float[nj]): nonuniform source point x-coords, in R
            y     (float[nj]): nonuniform source point y-coords, in R
            c     (complex[nj]): source strengths
            isign (int): if >=0, uses + sign in exponential, otherwise - sign
            eps   (float): precision requested (>1e-16)
            s     (float[nk]): nonuniform target x-frequencies, in R
            t     (float[nk]): nonuniform target y-frequencies, in R
            f     (complex[nk]): output values at target frequencies.

        """
        vis = np.empty((self.size_nu,), dtype=np.complex128, order='F')
        img = np.array(img, dtype=np.complex128, order='F')
        code = finufftpy.nufft2d3(self.x, self.y, img, isign=-1, eps=1e-6,
                                  s=self.u, t=self.v, f=vis, upsampfac=1.25,
                                  debug=0, spread_debug=0, spread_sort=2,
                                  fftw=1)
        if code == 1:
            raise Exception("eps too small in FT")
        return vis
