
* **Algorithm Options:** 
 * Method: ``(type: list)`` Short description
 * NL process: ``(type: list)`` Short description
 * Alpha (rad): ``(type: float)`` Short description
 * Gamma (Hz): ``(type: float)`` Short description
 * Material: ``(type: list)`` Short description
 * Dscan Parameter Scan: ``(type: group)`` Short description
  * Insertion min (mm): ``(type: float)`` Short description
  * Insertion max (mm): ``(type: float)`` Short description
  * Insertion step (mm): ``(type: float)`` Short description
 * MIIPS Parameter Scan: ``(type: group)`` Short description
  * Phase min (rad): ``(type: float)`` Short description
  * Phase max (rad): ``(type: float)`` Short description
  * Phase setp (rad): ``(type: float)`` Short description
* **Data Info** 
 * Loaded file: ``(type: text)`` Short description
 * Loaded node: ``(type: str)`` Short description
 * Trace Info ``(type: group)`` Short description
  * Wl0 (nm) ``(type: float)`` Short description
  * FWHM (nm) ``(type: float)`` Short description
  * Param Size ``(type: int)`` Short description
  * Wavelength Size ``(type: int)`` Short description
  * Wavelength scaling ``(type: float)`` Short description
  * Parameter scaling ``(type: float)`` Short description
 * Spectrum Info ``(type: group)`` Short description
  * Wl0 (nm) ``(type: float)`` Short description
  * FWHM (nm) ``(type: float)`` Short description
  * Fourier transform duration (fs) ``(type: float)`` Short description
  * Wavelength Size ``(type: int)`` Short description
  * Wavelength scaling ``(type: float)`` Short description
* **Processing** 
 * Grid settings: ``(type: group)`` Short description
  * lambda0 (nm): ``(type: float)`` Short description
  * Npoints: ``(type: list)`` Short description
  * Time resolution (fs): ``(type: float)`` Short description
 * Trace limits: ``(type: group)`` Short description
  * Crop Trace?: ``(type: bool)`` Short description
  * x0: ``(type: int)`` Short description
  * y0: ``(type: int)`` Short description
  * width: ``(type: int)`` Short description
  * height: ``(type: int)`` Short description
 * Substract trace background: ``(type: group)`` Short description
  * Substract?: ``(type: bool)`` Short description
  * wl0: ``(type: float)`` Short description
  * wl1: ``(type: float)`` Short description
 * Substract spectrum background: ``(type: group)`` Short description
  * Substract?: ``(type: bool)`` Short description
  * wl0: ``(type: float)`` Short description
  * wl1: ``(type: float)`` Short description
 * Process Spectrum ``(type: action)`` Short description
 * Process trace ``(type: action)`` Short description
 * Process Both ``(type: action)`` Short description
* **Retrieving** 
 * Algo type: ``(type: list)`` Short description
 * Verbose Info: ``(type: bool)`` Short description
 * Max iteration: ``(type: int)`` Short description
 * Uniform spectral response: ``(type: bool)`` Short description
 * Keep spectral intensity fixed ``(type: bool)`` Short description
 * Initial guess: ``(type: list)`` Short description
 * Initial Pulse Guess ``(type: group)`` Short description
  * FWHM (fs): ``(type: float)`` Short description
  * Phase amp. (rad): ``(type: float)`` Short description
 * Start Retrieval ``(type: action)`` Short description
 * Stop Retrieval ``(type: action)`` Short description
 * Propagate result ``(type: action)`` Short description