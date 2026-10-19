# -*- coding: utf-8 -*-

from sympy import symbols, Integral, exp, log, Add
from sympy.integrals.meijerint import meijerint_definite, meijerint_indefinite

def _make_integral_01():
    # Scaled down version of an integration by parts example
    params = symbols('u0 k l t'.split(), positive=True)
    u0, k, l, t = params
    M = exp(l*t)
    Q = k*u0*exp(-k*t)
    integral = Integral(M*Q, (t, 0, t))
    return integral, params


class TimeIntegration01:

    def setup(self):
        self.integral, self.params = _make_integral_01()
        self.values = {}
        self._set_ref()

    def _set_ref(self):
        u0, k, l, t = self.params
        self.ref = {
            'time_meijerint_definite': -k*u0/(k-l)*(exp(t*(l-k)) - 1),
        }

    def teardown(self):
        f, (t, _, _) = self.integral.function, self.integral.limits[0]
        for key, val in self.values.items():
            if key == 'time_meijerint_indefinite':
                # An antiderivative is only determined up to a constant
                if (val.diff(t) - f).simplify() != 0:
                    raise ValueError("Incorrect antiderivative, invalid"
                                     " timing: %s" % val)
                continue
            try:
                ref = self.ref[key]
            except KeyError:
                pass  # don't fail because of non-existent ref.
            else:
                if (ref - val).simplify() != 0:
                    raise ValueError("Incorrect result, invalid timing:"
                                     " %s != %s" % (ref, val))

    def time_doit(self):
        self.values['time_doit'] = self.integral.doit()

    def time_meijerint_definite(self):
        f, (t, a, b) = self.integral.function, self.integral.limits[0]
        # The convergence condition is not checked
        self.values['time_meijerint_definite'] = meijerint_definite(f, t, a, b)[0]

    def time_meijerint_indefinite(self):
        f, (t, _, _) = self.integral.function, self.integral.limits[0]
        self.values['time_meijerint_indefinite'] = meijerint_indefinite(f, t)


class _TimeIntegrationRisch(object):
    params = [1, 10, 100]

    def setup(self, n):
        x = symbols('x')
        self.integral = Integral(self.make_expr(n), x)
        self.values = {}

    def teardown(self, n):
        x = symbols('x')
        for key, val in self.values.items():
            ref = self.make_expr(n)
            if (ref - val.diff(x)).simplify() != 0:
                raise ValueError("Incorrect result, invalid timing:"
                                 " %s != %s" % (ref, val))

    def time_doit(self, n):
        self.values['time_doit'] = self.integral.doit()

    def time_doit_risch(self, n):
        self.values['time_doit_risch'] = self.integral.doit(risch=True)


class TimeIntegrationRisch01(_TimeIntegrationRisch):
    def make_expr(self, n):
        x = symbols('x')
        expr = x**n*exp(x)
        return expr

class TimeIntegrationRisch02(_TimeIntegrationRisch):
    def make_expr(self, n):
        x = symbols('x')
        expr = log(Add(*[exp(i*x) for i in range(n)])).diff(x)
        return expr

class TimeIntegrationRisch03(_TimeIntegrationRisch):
    def make_expr(self, n):
        x = symbols('x')
        expr = Add(*[log(x)**i for i in range(n)])
        return expr
