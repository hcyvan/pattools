from .sun.sun import deconvolution_sun, SunMarkers
from .moss.moss import deconvolution_moss, MossMarkers
from .loyfer.loyfer import deconvolution_loyfer, LoyferMarkers


def print_cell_type_helper(alg=None):
    if alg == 'sun':
        marker = SunMarkers()
        print("The deconvolution algorithm from Sun et al")
    elif alg == 'moss':
        marker = MossMarkers()
        print("The deconvolution algorithm from Moss et al")
    elif alg == 'loyfer':
        marker = LoyferMarkers()
        print("The deconvolution algorithm from Loyfer et al")
    else:
        raise Exception('unknown deconvolution algorithm')
    print(f"Cell/Tissue\t\tDetail")
    info = marker.get_cell_type_info()
    for k, v in info.items():
        print(f"{k}\t\t{v}")


__all__ = ['deconvolution_sun', 'deconvolution_moss', 'deconvolution_loyfer', 'print_cell_type_helper']
