from tvb.datatypes.equations import TemporalApplicableEquation
from tvb.basic.neotraits.api import Attr, Final

class SinglePulse(TemporalApplicableEquation):
    """
    A pulse train , offset with respect to the time axis.

    **Parameters**:
    
    * dT     :  pulse repetition period
    * onset         :  time of the first pulse
    * amp           :  amplitude of the pulse
    unit: milliseconds
    """

    equation = Final(
        label="Pulse Train",
        default="where((var>=onset)&((var-onset) < dT), amp, 0)",
        doc=""":math:`\\left\\{{\\begin{array}{rl}amp,&{\\text{if }} ((var-onset) < dT \\and var > onset\\\\0,&{\\text{otherwise }}\\end{array}}\\right.`""")

    parameters = Attr(
        field_type=dict,
        default=lambda: {"dT": 0.1, "amp": 0.001, "onset": 500},
        label="Pulse Train Parameters")
