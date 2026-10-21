"""
Reading, writing, plotting, cleaning and continuum normalization of 1D spectra
of (massive) stars.

Conventions shared by every function in this module:
    - A spectrum is a pair (wave, flux) of 1D float NumPy arrays with the
      same length, or a spectrum file (ASCII or FITS) read with 'readspec'.
    - Wavelengths are kept in the units of the input spectrum (usually Å)
      and are never converted or corrected (e.g. for the barycentric
      velocity; see 'panorama.obs.bary_corr'); plot axes are labeled in Å.
    - File arguments are paths (path + filename), given either as strings
      or as pathlib.Path objects. Output directories are never created.
    - Spectra are written with 'writespec' (ASCII or FITS, chosen from the
      extension of the output file), so every spectrum written by the module
      can be read back with 'readspec'.
    - All functions validate their inputs before computing.
    - Interactive functions require an interactive Matplotlib backend (e.g.
      QtAgg, TkAgg; '%matplotlib qt' in Jupyter).
"""

from __future__ import annotations

# =====================================================================
# ============================== SPECTRA ==============================
# =====================================================================

# ==============================
def readspec(
             file : str | Path
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Read a 1D spectrum file (ASCII or FITS, optionally gzip-compressed) and return its wavelength and flux arrays.

    ----------
    Parameters
    ----------
        file : str | Path
            Path and filename of the spectrum file. The format is detected
            from the extension (a final '.gz' is accepted).

    -------
    Returns
    -------
        wave : np.ndarray[float]
            Wavelength array (strictly increasing), in the units of the file.

        flux : np.ndarray[float]
            Flux array (same length as 'wave').

    ------
    Raises
    ------
        TypeError
            If 'file' is not a path.

        FileNotFoundError
            If the file does not exist.

        IsADirectoryError
            If 'file' is a directory.

        PermissionError
            If the file cannot be read.

        ValueError
            If the extension of the file is not supported, or the file does
            not contain a valid spectrum (no data, missing WCS keywords, arrays
            of different lengths, non-finite or non-monotonic wavelengths...).

        RuntimeError
            If the file cannot be read for any other reason (e.g. a corrupt
            FITS file).

    -----
    Notes
    -----
        1. ASCII FORMATS (.txt, .dat, .asc, .ascii)
            - Two or more columns (wavelength, flux, optional extra columns).
              Only the first two must be numeric; extra columns (e.g. quality
              flags) may contain text and trigger a warning.

            - Lines starting with '#' and empty lines are ignored.
              Non-numeric lines before the first data line (column names,
              e.g. 'wavelength flux') are skipped with a warning.

            - Separators can be whitespace, commas or semicolons.

            - Every data line must have the same number of columns.

        2. FITS FORMATS (.fits, .fit, .fts)
            - Binary table with wavelength and flux columns (commonly named
              WAVE / WAVELENGTH / LAMBDA / LAM and FLUX / INTENSITY / FLX /
              SPEC / FLAMBDA / COUNTS / CTS), searched across every extension.
              Tables without recognizable columns are skipped.

            - Otherwise, image with WCS keywords (CRVAL1, CDELT1 or CD1_1,
              CRPIX1) in the PrimaryHDU or, if that one is empty, in the first
              image extension with data. For multi-dimensional images (e.g.
              flux + error, IRAF multispec) the first row is used; a missing
              CRPIX1 means pixel 1 (FEROS / IRAF convention).

            - Logarithmic scales: DC-FLAG = 1 (IRAF, log10), CTYPE1 ending in
              '-LOG' (FITS standard; Greisen et al. 2006), or ln(lambda)
              values in files with 'log' in their name (e.g. HERMES
              '*_log_merged_cf.fits' spectra).

        3. CHECKS
            - Descending wavelength arrays are reversed. Non-finite fluxes
              are kept, with a warning.

        4. READERS
            - The format is chosen by '_spectrum_file_format' and each one is
              read by a private function of 'support/spec/spec_support.py'
              ('_read_ascii', '_read_fits'), which only returns the raw
              arrays: the generic checks above are done here for every
              format.

            - ValueError raised by a reader (invalid content) is propagated
              unchanged; any other exception is chained to a RuntimeError
              with the file name.

            - Every file written by 'writespec' can be read back with this
              function, since both use the same extensions.

    ----------
    References
    ----------
        - Pence, W. D., Chiappetti, L., Page, C. G. et al. 2010, A&A, 524, A42.
        - Greisen, E. W., Calabretta, M. R., Valdes, F. G. et al. 2006,
          A&A, 446, 747.
        - Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022,
          ApJ, 935, 167.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import numpy as np
    import warnings


    # ##### PANORAMA #####
    from .support.spec.spec_support import _read_ascii, _read_fits, _spectrum_file_format
    from .support.tex.tex_support import _check_input_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    file           = _check_input_file(file, "file", "Spectrum file", "spectrum file")
    file_format, _ = _spectrum_file_format(file, name="file", label="Spectrum file")  # 'ascii' or 'fits', from the extension (unsupported ones are rejected)


    # ####################  COMPUTATION  ####################
    # ##########  READING  ##########
    try:
        wave, flux = _read_ascii(file) if file_format == "ascii" else _read_fits(file)  # Raw arrays; the generic checks are done below
    except ValueError:
        raise  # Invalid content of the file (no data, missing WCS keywords...): propagated unchanged
    except Exception as e:
        raise RuntimeError(f"Spectrum file could not be read: '{file}' ({type(e).__name__}: {e}).") from e
    wave, flux = np.asarray(wave, dtype=float), np.asarray(flux, dtype=float)  # Same dtype whatever the reader


    # ##########  CHECKS  ##########
    if wave.ndim != 1 or flux.ndim != 1:
        raise ValueError(f"Wavelength and flux arrays of spectrum file '{file}' must both be 1D. Got shapes {wave.shape} and {flux.shape}.")
    elif wave.shape != flux.shape:
        raise ValueError(f"Wavelength and flux arrays of spectrum file '{file}' must have the same length. Got {wave.size} and {flux.size}.")
    elif wave.size == 0:
        raise ValueError(f"Spectrum file '{file}' contains no spectral data.")
    #
    wave_bad = np.count_nonzero(~np.isfinite(wave))  # A wavelength grid with NaN / inf cannot be used
    if wave_bad:
        raise ValueError(f"Wavelength array of spectrum file '{file}' contains {wave_bad} non-finite value(s).")
    #
    flux_bad = np.count_nonzero(~np.isfinite(flux))  # Non-finite fluxes (e.g. masked pixels) are kept
    if flux_bad:
        warnings.warn(f"Flux array of spectrum file '{file}' contains {flux_bad} non-finite value(s).", stacklevel=2)


    # ##########  ORDERING  ##########
    if wave.size > 1 and np.all(np.diff(wave) < 0):  # Descending scale (e.g. negative CDELT1)
        wave, flux = wave[::-1].copy(), flux[::-1].copy()  # Copies: contiguous arrays instead of reversed views
    if np.any(np.diff(wave) <= 0):  # Repeated or unsorted wavelengths after the reversal
        raise ValueError(f"Wavelength array of spectrum file '{file}' is not strictly increasing (overlapping echelle orders or duplicated points?).")


    # ####################  OUTPUT  ####################
    return wave, flux




# ==============================
def writespec(
              spec        : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path,
              output_file : str | Path,
              overwrite   : bool = False,
              header      : str | None = None
    ) -> Path:
    """
    Write a 1D spectrum to an ASCII or FITS file (format chosen from the extension of 'output_file') and return its path.

    ----------
    Parameters
    ----------
        spec : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path
            Spectrum to write: a (wave, flux) tuple or [wave, flux] list of 1D
            array-like objects with the same length, or the path to a
            spectrum file readable by 'readspec' (format conversion).

        output_file : str | Path
            Path and filename of the output spectrum file. The format is
            chosen from its extension (case-insensitive): ASCII (.txt, .dat,
            .asc, .ascii) or FITS (.fits, .fit, .fts); a final '.gz' writes it
            gzip-compressed (e.g. 'HD886.fits.gz'). Its directory must exist.

        overwrite : bool, optional
            If True, an existing 'output_file' is replaced. Default is False.

        header : str | None, optional
            Text stored with the spectrum, e.g. a description of how it was
            processed: one '#' comment line per line of the text in ASCII
            files (give it without the '#' prefix), or one COMMENT card per
            line in FITS files (ASCII characters only). Default is None (no
            extra lines).

    -------
    Returns
    -------
        output_file : Path
            Path of the written spectrum file.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the spectrum is not valid (arrays not 1D, of different lengths,
            empty, with non-finite or non-monotonic wavelengths), the
            extension of 'output_file' is not supported, or 'header' contains
            non-ASCII characters for a FITS file.

        IsADirectoryError
            If 'output_file' is an existing directory.

        FileNotFoundError
            If the directory of 'output_file' does not exist.

        NotADirectoryError
            If the parent of 'output_file' is not a directory.

        FileExistsError
            If 'output_file' already exists and 'overwrite' is False.

        RuntimeError
            If the file cannot be written for any other reason (e.g. no write
            permission).

    -----
    Notes
    -----
        1. ASCII FORMAT (.txt, .dat, .asc, .ascii)
            - Two tab-separated columns (wavelength, flux) after the '#'
              header lines: the lines of 'header' (if given) and the column
              names 'wave' and 'flux' (also tab-separated). Wavelengths are
              written with '%.12g' and fluxes with '%.10g'; non-finite fluxes
              are written as 'nan' / 'inf'.

        2. FITS FORMAT (.fits, .fit, .fts)
            - Empty PrimaryHDU plus a binary table extension 'SPECTRUM' with
              the float64 columns WAVE and FLUX. No units are written
              (TUNITn), because wavelengths keep the units of the input.

            - The lines of 'header' (if given) are stored as COMMENT cards of
              the 'SPECTRUM' extension; FITS headers only accept ASCII
              characters.

        3. SPECTRUM
            - Descending wavelength arrays are reversed, so that the file is
              readable by 'readspec'. Repeated or unsorted wavelengths are
              rejected. Non-finite fluxes are kept, with a warning.

            - The extensions are those of 'readspec' ('_ASCII_EXTENSIONS' and
              '_FITS_EXTENSIONS', see 'support/spec/spec_support.py'), so
              every written file can be read back with 'readspec'.

            - Spectrum files raise the exceptions of 'readspec'.

        4. OTHER FUNCTIONS
            - 'rmvspec' and 'normspec' write their output spectra with this
              function, replacing the existing file (they save repeatedly
              from their interactive windows); 'normspec' stores the
              description of the continuum fit in 'header'.

            - Example (FITS to gzip-compressed ASCII):
                >>> writespec("HD886.fits", "HD886.txt.gz")

    ----------
    References
    ----------
        - Pence, W. D., Chiappetti, L., Page, C. G. et al. 2010, A&A, 524, A42.
        - Astropy Collaboration, Price-Whelan, A. M., Lim, P. L. et al. 2022,
          ApJ, 935, 167.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    from astropy.io import fits
    import numpy as np
    import warnings


    # ##### PANORAMA #####
    from .support.spec.spec_support import _load_spectrum, _spectrum_file_format
    from .support.tex.tex_support import _check_output_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    wave, flux = _load_spectrum(spec)  # (wave, flux) pair or spectrum file (format conversion)
    #
    output_file    = _check_output_file(output_file, allow_none=False)
    file_format, _ = _spectrum_file_format(output_file)  # 'ascii' or 'fits', from the extension
    #
    if not isinstance(overwrite, bool):
        raise TypeError(f"Overwrite flag 'overwrite' must be a boolean. Got {type(overwrite).__name__}.")
    #
    if header is not None:
        if not isinstance(header, str):
            raise TypeError(f"Spectrum file header 'header' must be a string or None. Got {type(header).__name__}.")
        elif file_format == "fits" and not header.isascii():
            raise ValueError(f"Spectrum file header 'header' must contain only ASCII characters in a FITS file. Got '{header}'.")


    # ##########  EXISTING FILE  ##########
    if output_file.exists() and not overwrite:
        raise FileExistsError(f"Output file already exists (use 'overwrite=True' to replace it): '{output_file}'.")


    # ##########  HEADER LINES  ##########
    header_lines = header.splitlines() if header is not None else []  # One '#' line (ASCII) or COMMENT card (FITS) per line


    # ####################  COMPUTATION  ####################
    # ##########  CHECKS  ##########
    wave_bad = np.count_nonzero(~np.isfinite(wave))  # The file could not be read back by 'readspec'
    if wave_bad:
        raise ValueError(f"Wavelength array of the spectrum 'spec' must be finite. Got {wave_bad} non-finite value(s).")
    #
    flux_bad = np.count_nonzero(~np.isfinite(flux))
    if flux_bad:
        warnings.warn(f"Flux array of the spectrum 'spec' contains {flux_bad} non-finite value(s); they are written as they are.", stacklevel=2)


    # ##########  ORDERING  ##########
    if wave.size > 1 and np.all(np.diff(wave) < 0):  # Descending scale: reversed so that 'readspec' accepts the file
        wave, flux = wave[::-1].copy(), flux[::-1].copy()
    wave_steps = np.count_nonzero(np.diff(wave) <= 0)  # Repeated or unsorted wavelengths, rejected by 'readspec'
    if wave_steps:
        raise ValueError(f"Wavelength array of the spectrum 'spec' must be strictly monotonic. Got {wave_steps} repeated or unsorted wavelength(s).")


    # ####################  OUTPUT  ####################
    # ##########  WRITE OUTPUT  ##########
    try:


        # ##### OPTION 1: ASCII #####
        if file_format == "ascii":
            np.savetxt(
                       output_file,  # A '.gz' name is compressed by NumPy
                       np.column_stack([wave, flux]),  # Two columns: wavelength, flux
                       fmt=["%.12g", "%.10g"],  # Enough digits for high-resolution wavelengths and absolute fluxes
                       delimiter="\t",
                       header="\n".join(header_lines + ["wave\tflux"]),  # Column names in the last '#' line
                       comments="# ",
                       encoding="utf-8"
                      )


        # ##### OPTION 2: FITS #####
        else:
            table_hdu = fits.BinTableHDU.from_columns(
                                                      [
                                                       fits.Column(name="WAVE", format="D", array=wave),  # 'D' = float64
                                                       fits.Column(name="FLUX", format="D", array=flux)
                                                      ],
                                                      name="SPECTRUM"
                                                     )
            for line in header_lines:
                table_hdu.header.add_comment(line)  # One COMMENT card per line of 'header'
            fits.HDUList([fits.PrimaryHDU(), table_hdu]).writeto(output_file, overwrite=overwrite)  # Empty primary HDU + 'SPECTRUM' table; '.gz' names are compressed by astropy
    except OSError as e:  # E.g. no write permission
        raise RuntimeError(f"Spectrum file could not be written: '{output_file}'.") from e


    # ##########  RESULTS  ##########
    return output_file




# ==============================
def plotspec(
             spec       : list[tuple[np.ndarray[float], np.ndarray[float]] | str | Path] | tuple[np.ndarray[float], np.ndarray[float]] | str | Path,
             wave_lim   : tuple[float, float] | None = None,
             show_lines : bool = True,
             show       : bool = True
    ) -> Figure:
    """
    Plot one or several spectra in an interactive window with check buttons to toggle their visibility and the main O / B star lines.

    ----------
    Parameters
    ----------
        spec : list[tuple[np.ndarray[float], np.ndarray[float]] | str | Path] | tuple[np.ndarray[float], np.ndarray[float]] | str | Path
            Spectra to plot. A single spectrum is a (wave, flux) tuple of 1D
            array-like objects with the same length, or the path to a spectrum
            file readable by 'readspec'. Several spectra are given as a list
            of such items.

        wave_lim : tuple[float, float] | None, optional
            Wavelength range (min, max) shown, in the units of the spectra.
            The flux axis is scaled to the points inside it. Default is None
            (automatic limits).

        show_lines : bool, optional
            If True, the main diagnostic lines of O and B stars (see Notes)
            inside the wavelength range of the spectra are drawn as vertical
            dotted lines, labeled at the bottom of the axes. They can be
            hidden / shown with the 'Lines (OB-type stars)' check button.
            Default is True.

        show : bool, optional
            If True, display the figure with plt.show(). Default is True.

    -------
    Returns
    -------
        fig : Figure
            The matplotlib.figure.Figure with the spectra and the check
            buttons (not saved; see 'savefig_pkl').

    ------
    Raises
    ------
        TypeError
            If 'spec' or one of its items has the wrong type, 'wave_lim' is
            not a pair of real numbers, or 'show_lines' or 'show' is not a
            boolean.

        ValueError
            If 'spec' is an empty list, a spectrum is not valid or 'wave_lim'
            is not a valid range.

    -----
    Notes
    -----
        1. SPECTRA
            - Spectrum files are labeled with their file name (without
              extension) and (wave, flux) tuples as 'Spectrum N'. Duplicated
              labels get a ' (k)' suffix.

            - The wavelength axis is labeled in Å, the usual unit of the
              spectra; the data are not converted.

            - With a non-interactive backend (e.g. Agg) 'plt.show' does
              nothing and the check buttons do not respond; with the inline
              backend of Jupyter the figure is shown in the notebook.

            - The check buttons keep working while the returned figure is
              open: a reference to the widget is stored in the figure.

            - Spectrum files raise the exceptions of 'readspec'.

        2. SPECTRAL LINES
            - 'support/spec/spec_OBlines.py' contains 40 diagnostic lines of O
              and B stars: Balmer lines, He I, He II, and C, N, O, Mg and Si
              lines used in spectral classification. Wavelengths are
              rest-frame air wavelengths in Å.

            - Each label is drawn in a blended transform: its x coordinate is
              in wavelength units and its y coordinate in axes units, so it
              follows its line when zooming / panning in wavelength and always
              stays at the bottom of the axes when zooming in flux. Labels of
              lines outside the visible range are clipped.

            - Colors: H I red, He I blue, He II purple; metals by element
              (C brown, N orange, O cyan, Mg olive, Si green).

            - The lines are drawn at rest.

            - Blends are labeled by their main component (e.g. 'He I+II 4026',
              'C III 4650' for C III 4647-4651).

    ----------
    References
    ----------
        - Walborn, N. R. & Fitzpatrick, E. L. 1990, PASP, 102, 379.
        - Gray, R. O. & Corbally, C. J. 2009,
          Stellar Spectral Classification (Princeton University Press).
        - Sota, A., Maíz Apellániz, J., Walborn, N. R. et al. 2011,
          ApJS, 193, 24.
        - Hunter, J. D. 2007, CSE, 9, 90.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PYTHON #####
    import matplotlib.pyplot as plt
    from matplotlib.widgets import CheckButtons
    import numpy as np
    import os
    from pathlib import Path


    # ##### PANORAMA #####
    from .support.spec.spec_support import _check_wave_lim, _draw_spectral_lines, _load_spectrum, _SpectrumToggle


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    if not isinstance(spec, (list, tuple, str, os.PathLike)):
        raise TypeError(f"Spectra to plot 'spec' must be a (wave, flux) tuple, a path (str or pathlib.Path) or a list of them. Got {type(spec).__name__}.")
    items = spec if isinstance(spec, list) else [spec]  # A single spectrum is treated as a list of one
    if not items:
        raise ValueError("Spectra to plot 'spec' must contain at least one spectrum. Got an empty list.")
    for i, item in enumerate(items):
        if not isinstance(item, (tuple, str, os.PathLike)):
            raise TypeError(f"Spectrum 'spec[{i}]' must be a (wave, flux) tuple or a path (str or pathlib.Path) to a spectrum file. Got {type(item).__name__}.")
    #
    _check_wave_lim(wave_lim)  # None or a finite (min, max) pair with min < max
    #
    if not isinstance(show_lines, bool):
        raise TypeError(f"Spectral lines flag 'show_lines' must be a boolean. Got {type(show_lines).__name__}.")
    #
    if not isinstance(show, bool):
        raise TypeError(f"Display flag 'show' must be a boolean. Got {type(show).__name__}.")


    # ####################  PARAMETERS  ####################
    lines_label = "Lines (OB-type stars)"  # Check button of the spectral lines (reserved label)


    # ####################  COMPUTATION  ####################
    # ##########  READING  ##########
    spectra = []
    labels  = []
    for i, item in enumerate(items):
        wave, flux = _load_spectrum(item, name=f"spec[{i}]")  # The index appears in the error messages
        label      = Path(item).stem if isinstance(item, (str, os.PathLike)) else f"Spectrum {i + 1}"
        #
        base_label, counter = label, 1
        while label in labels or label == lines_label:  # Unique labels: the check buttons identify spectra by label
            label    = f"{base_label} ({counter})"
            counter += 1  # ' (1)', ' (2)'... until the label is free
        spectra.append((wave, flux))
        labels.append(label)


    # ##########  WAVELENGTH RANGE  ##########
    wave_min = min(np.nanmin(wave) for wave, _ in spectra)  # Range where the spectral lines are drawn
    wave_max = max(np.nanmax(wave) for wave, _ in spectra)


    # ####################  PLOTTING  ####################
    # ##########  FIGURE  ##########
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(left=0.08, right=0.75, top=0.9, bottom=0.12)  # Right margin left free for the check buttons


    # ##########  SPECTRA  ##########
    spectra_plots = [ax.plot(wave, flux, label=label, linewidth=1.2, alpha=0.85)[0] for (wave, flux), label in zip(spectra, labels)]  # One Line2D per spectrum, labeled as its check button


    # ##########  SPECTRAL LINES  ##########
    line_artists = _draw_spectral_lines(ax, wave_min, wave_max) if show_lines else []  # Two artists per line (vertical line + label)


    # ##########  AXES  ##########
    if wave_lim is not None:
        ax.set_xlim(wave_lim)
        flux_in = [flux[(wave >= wave_lim[0]) & (wave <= wave_lim[1])] for wave, flux in spectra]  # Fluxes inside the range, per spectrum
        flux_in = np.concatenate(flux_in)
        flux_in = flux_in[np.isfinite(flux_in)]  # NaN / inf would break the limits
        if flux_in.size:  # Otherwise the flux axis would still span the whole spectra
            flux_min, flux_max = flux_in.min(), flux_in.max()
            pad                = 0.05 * (flux_max - flux_min) if flux_max > flux_min else (0.05 * abs(flux_max) or 1.)  # 5% margin; constant or zero flux also handled
            ax.set_ylim(flux_min - pad, flux_max + pad)
    ax.set_xlabel("Wavelength [Å]", fontsize=12, fontweight="normal")
    ax.set_ylabel("Flux", fontsize=12, fontweight="normal")
    ax.tick_params(axis="both", which="major", labelsize=10)
    ax.minorticks_on()
    ax.grid(True, which="major", linestyle="-", alpha=0.5)  # Solid major grid, dotted minor grid
    ax.grid(True, which="minor", linestyle=":", alpha=0.2)


    # ##########  CHECK BUTTONS  ##########
    check_labels  = labels + ([lines_label] if line_artists else [])  # Spectral lines button only if any line was drawn
    widget_height = min(0.06 * len(check_labels), 0.8)  # Figure fraction: grows with the number of buttons, up to 80%
    rax           = fig.add_axes([0.76, 0.5 - widget_height / 2, 0.2, widget_height])  # [left, bottom, width, height], vertically centered
    for spine in rax.spines.values():
        spine.set_visible(False)  # Buttons without a frame
    #
    check               = CheckButtons(rax, check_labels, [True] * len(check_labels))  # Every item visible at start
    toggle              = _SpectrumToggle(fig, spectra_plots, line_artists, lines_label)
    fig._plotspec_check = (check, toggle)  # Keep references: an unreferenced widget or callback stops responding
    for label_text, line in zip(check.labels, spectra_plots):
        label_text.set_color(line.get_color())  # Each button label takes the color of its spectrum
        label_text.set_fontweight("bold")
    if line_artists:
        check.labels[-1].set_color("0.35")  # Gray label for the spectral lines button (last one)
    check.on_clicked(toggle._on_clicked)


    # ##########  SHOW PLOT  ##########
    if show:
        plt.show()  # Inline backend: displays the figure; interactive backends: opens a window


    # ####################  OUTPUT  ####################
    return fig




# ==============================
def rmvspec(
            spec        : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path,
            output_file : str | Path | None = None,
            wave_lim    : tuple[float, float] | None = None,
            autosave    : bool = True,
            show_lines  : bool = True
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Remove regions of a spectrum interactively by mouse selection and return the cleaned spectrum.

    ----------
    Parameters
    ----------
        spec : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path
            Spectrum to clean: a (wave, flux) tuple or [wave, flux] list of 1D
            array-like objects with the same length, or the path to a
            spectrum file readable by 'readspec'.

        output_file : str | Path | None, optional
            File where the cleaned spectrum is written with 'writespec' (ASCII
            or FITS, chosen from its extension; see 'writespec'), with 's' /
            'S' and, if 'autosave', when the window is closed. Exported spans
            are written next to it as '<stem>_rmvspans.txt' /
            '<stem>_selectedspans.txt' and the figure as '<stem>.png' (with a
            '<stem>.pkl' copy, see 'savefig_pkl'), where '<stem>' is
            'output_file' without its extension (and '.gz'). Its directory
            must exist. If None, nothing is written and the cleaned spectrum
            is only returned. Default is None.

        wave_lim : tuple[float, float] | None, optional
            Wavelength range (min, max) to display and clean. Points outside
            it are excluded from the tool and from the returned spectrum.
            Default is None (whole spectrum).

        autosave : bool, optional
            If 'output_file' is given, write the cleaned spectrum when the
            window is closed, unless it has not changed since the last save.
            Default is True.

        show_lines : bool, optional
            If True, the main diagnostic lines of O and B stars ('_OB_LINES',
            from 'spec_OBlines') inside the spectrum are drawn as vertical
            dotted lines, labeled at the bottom of the axes. They can be
            hidden / shown with 'l'. Default is True.

    -------
    Returns
    -------
        wave_clean : np.ndarray[float]
            Wavelengths of the kept points.

        flux_clean : np.ndarray[float]
            Fluxes of the kept points.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If the spectrum or 'wave_lim' are not valid, fewer than two finite
            points remain inside 'wave_lim', or the extension of
            'output_file' is not supported by 'writespec'.

        IsADirectoryError
            If 'output_file' is an existing directory.

        FileNotFoundError
            If the directory of 'output_file' does not exist.

        NotADirectoryError
            If the parent of 'output_file' is not a directory.

    -----
    Notes
    -----
        1. KEYBOARD SHORTCUTS (they take priority over Matplotlib's)
            - r : remove all selected (red) spans.

            - x : discard the last selected span.

            - c : clear all selected spans.

            - a : restore all removed points.

            - u / U : undo the last edit / redo.

            - g : toggle the grid.

            - l : show / hide the O / B star lines (it replaces Matplotlib's
              log-scale toggle in this window).

            - h : print the keyboard help.

            - s / S : save the cleaned spectrum to 'output_file' / and a PNG
              of the figure (plus its '.pkl' copy).

            - d : export the removed spans (merged) to '<stem>_rmvspans.txt'.

            - D : export the pending selected spans to
              '<stem>_selectedspans.txt'.

        2. MATPLOTLIB
            - Matplotlib shortcuts not used by the tool still work: 'o'
              (zoom), 'p' (pan), 'q' (close), arrows (back / forward)... The
              zoom is kept after every removal. No global Matplotlib setting
              is modified.

            - An interactive backend is required (e.g. QtAgg, TkAgg; in
              Jupyter, '%matplotlib qt'). With a non-interactive one (e.g.
              Agg, the inline backend) a warning is issued, the figure is
              closed without autosaving and the spectrum is returned
              unmodified (within 'wave_lim').

        3. INPUT SPECTRUM
            - Non-finite points are discarded (with a warning) and the
              spectrum is sorted by wavelength if needed.

        4. OUTPUT FILES
            - The cleaned spectrum is written with 'writespec', replacing the
              previous file, so it can be read back with 'readspec'.

            - Errors while saving inside the window (e.g. permissions, or a
              spectrum with repeated wavelengths, which 'writespec' rejects)
              are printed, not raised.

        5. SPECTRAL LINES
            - The lines are the same as in 'plotspec': rest-frame air
              wavelengths in Å from 'spec_OBlines', colored by species. Their
              labels always stay at the bottom of the axes when zooming.

            - Spectra not corrected for the radial velocity of the star show
              their lines shifted by v_rad / c * λ from the marks.

            - Spectrum files raise the exceptions of 'readspec'.

    ----------
    References
    ----------
        - Hunter, J. D. 2007, CSE, 9, 90.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from .support.spec.spec_support import _check_wave_lim, _InteractiveSpectrumCleaner, _load_spectrum, _prepare_spectrum, _spectrum_file_format
    from .support.tex.tex_support import _check_output_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    wave, flux = _load_spectrum(spec)  # (wave, flux) pair or spectrum file
    #
    output_file = _check_output_file(output_file, allow_none=True)  # None: nothing is written
    if output_file is not None:
        _spectrum_file_format(output_file)  # Unsupported extensions are rejected before the window opens
    #
    _check_wave_lim(wave_lim)
    #
    if not isinstance(autosave, bool):
        raise TypeError(f"Autosave flag 'autosave' must be a boolean. Got {type(autosave).__name__}.")
    #
    if not isinstance(show_lines, bool):
        raise TypeError(f"Spectral lines flag 'show_lines' must be a boolean. Got {type(show_lines).__name__}.")


    # ##########  SPECTRUM  ##########
    wave, flux = _prepare_spectrum(wave, flux, wave_lim, keep_flux_nonfin=False)  # Non-finite fluxes can be neither drawn nor selected


    # ####################  COMPUTATION  ####################
    cleaner                = _InteractiveSpectrumCleaner(wave, flux, output_file, show_lines=show_lines, autosave=autosave)  # Builds the figure and takes the keyboard
    wave_clean, flux_clean = cleaner._run()  # Blocks until the window is closed


    # ####################  OUTPUT  ####################
    return wave_clean, flux_clean




# ==============================
def normspec(
             spec        : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path,
             output_file : str | Path | None = None,
             wave_lim    : tuple[float, float] | None = None,
             function    : str = "spline3",
             order       : int = 1,
             sample      : str | list[tuple[float, float]] | tuple[float, float] | None = "*",
             naverage    : int = 1,
             low_reject  : float = 2.,
             high_reject : float = 0.,
             niterate    : int = 10,
             grow        : float = 1.,
             output_type : str = "ratio",
             replace     : bool = False,
             markrej     : bool = True,
             interactive : bool = True
    ) -> tuple[np.ndarray[float], np.ndarray[float] | None]:
    """
    Fit the continuum of a 1D spectrum and normalize it, either automatically or interactively.

    ----------
    Parameters
    ----------
        spec : tuple[np.ndarray[float], np.ndarray[float]] | list[np.ndarray[float]] | str | Path
            Spectrum to normalize: a (wave, flux) tuple or [wave, flux] list
            of 1D array-like objects with the same length, or the path to a
            spectrum file readable by 'readspec'.

        output_file : str | Path | None, optional
            File where the output spectrum (wavelength, output flux) is
            written with 'writespec', in ASCII or FITS format (chosen from its
            extension; see 'writespec'), with a header describing the fit.
            Its directory must exist. In interactive mode it is written when
            'w' is pressed; in non-interactive mode, automatically. If None,
            nothing is written. Default is None.

        wave_lim : tuple[float, float] | None, optional
            Wavelength range (min, max) to normalize. Points outside it are
            excluded from the fit and from the returned spectrum. Default is
            None (whole spectrum).

        function : str, optional
            Continuum function: 'spline3' (cubic spline), 'spline1' (linear
            spline), 'legendre' or 'chebyshev' (polynomials). Default is
            "spline3".

        order : int, optional
            Number of spline pieces, or number of polynomial terms (degree +
            1). Default is 1.

        sample : str | list[tuple[float, float]] | tuple[float, float] | None, optional
            Wavelength ranges used in the fit, in the units of the spectrum:
            '*' or None (whole spectrum), an IRAF range string (e.g.
            '4000:4090,4150:4300'), a list of (min, max) tuples or a single
            (min, max) tuple. Default is "*".

        naverage : int, optional
            Number of consecutive pixels combined into each sample point:
            their average if positive, their median if negative. Non-zero.
            Default is 1.

        low_reject : float, optional
            Rejection limit below the fit, in units of the rms of the
            residuals. Zero disables it. Default is 2.

        high_reject : float, optional
            Rejection limit above the fit, in units of the rms of the
            residuals. Zero disables it. The defaults (2, 0) suit absorption
            line spectra. Default is 0.

        niterate : int, optional
            Maximum number of rejection iterations (>= 1). Default is 10.

        grow : float, optional
            Rejection growing radius, in sample points (pixels if naverage =
            1): points within this distance of a rejected point are also
            rejected. Only its integer part is used. Default is 1.

        output_type : str, optional
            Output spectrum: 'ratio' (flux / continuum), 'difference'
            (flux - continuum), 'fit' (continuum) or 'data' (flux). Unique
            abbreviations are accepted (e.g. 'r', 'di'). Default is "ratio".

        replace : bool, optional
            Replace the rejected and deleted pixels by the fit before
            computing the 'ratio', 'difference' and 'data' outputs. Default is
            False.

        markrej : bool, optional
            Mark the rejected points in the interactive plot. Default is True.

        interactive : bool, optional
            Open the interactive window (True) or just fit, save and return
            (False, useful for batch processing). Default is True.

    -------
    Returns
    -------
        wave : np.ndarray[float]
            Wavelength array (within 'wave_lim', if given).

        result : np.ndarray[float] | None
            Output spectrum ('output_type') of the last continuum fitted
            (automatic or manual), or None if no continuum could be fitted in
            interactive mode.

    ------
    Raises
    ------
        TypeError
            If any input has the wrong type.

        ValueError
            If any input has an invalid value (unknown 'function' or
            'output_type', non-positive 'order', zero 'naverage', negative
            rejection limits or 'grow', malformed 'sample', extension of
            'output_file' not supported by 'writespec'...), fewer than two
            points remain inside 'wave_lim', or, in non-interactive mode, the
            continuum cannot be fitted or written (e.g. repeated
            wavelengths).

        IsADirectoryError
            If 'output_file' is an existing directory.

        FileNotFoundError
            If the directory of 'output_file' does not exist.

        NotADirectoryError
            If the parent of 'output_file' is not a directory.

        RuntimeError
            If, in non-interactive mode, 'output_file' cannot be written for
            any other reason (e.g. no write permission).

    -----
    Notes
    -----
        1. AUTOMATIC FIT (IRAF 'continuum' / 'icfit')
            - Sample points are taken from the wavelength ranges in 'sample',
              combining 'naverage' consecutive pixels into each point.

            - The function is least-squares fitted to the sample points.

            - Points with residuals below -low_reject * rms or above
              high_reject * rms are rejected, together with the points within
              'grow' of them, and the function is refitted. This is repeated
              up to 'niterate' times or until no new points are rejected.
              Rejected points stay rejected.

            - The output spectrum is computed from the final continuum.
              Ratios where the continuum is <= 0 are set to NaN.

        2. INTERACTIVE CONTROLS: AUTOMATIC CONTINUUM
            - 'f': fit the continuum with the current parameters.

            - '+' / '-': increase / decrease the order.

            - '1' / '3': linear / cubic spline; 'l' / 'c': Legendre /
              Chebyshev polynomial.

            - 's' + 's': mark a sample region (press 's' at both ends); the
              first region marked replaces the whole-spectrum sample.

            - 't': reset the sample to the whole spectrum.

            - 'd' / 'u': delete / undelete the sample point nearest to the
              cursor.

        3. INTERACTIVE CONTROLS: MANUAL CONTINUUM AND GENERAL
            - Left-click: add a continuum point. Right-click on a point:
              remove it. Enter: fit a cubic spline through the manual points.

            - 'r': reset to the initial parameters; 'w': save the output
              spectrum to 'output_file'; 'q': close; '?': print the controls.

            - Matplotlib's default keyboard shortcuts are disabled in this
              figure (the toolbar buttons still work).

            - With a non-interactive backend (e.g. Agg, the inline backend) a
              warning is issued and the result of the initial automatic fit is
              returned (nothing is written); use 'interactive=False' for batch
              processing. With the inline backend of Jupyter, the figure of
              the automatic fit is shown in the notebook.

            - With an interactive backend the function returns when the window
              is closed, also in Jupyter / IPython (e.g. '%matplotlib qt'),
              where a plain 'plt.show' would return at once.

        4. INPUT SPECTRUM
            - Points with non-finite wavelengths are discarded (with a
              warning) and the spectrum is sorted by wavelength if needed.
              Non-finite fluxes are kept but never used as sample points.

        5. OUTPUT FILE
            - The output spectrum is written with 'writespec', replacing the
              previous file, so it can be read back with 'readspec'. The fit
              is described in '#' lines before the 'wave' / 'flux' column
              names (ASCII) or in COMMENT cards (FITS).

            - The status line and the header of 'output_file' always describe
              the continuum shown, even if a later fit failed.

            - Non-finite output fluxes (e.g. ratios where the continuum is
              <= 0) are written as they are, with a warning from 'writespec'.

            - Spectrum files raise the exceptions of 'readspec'.

    ----------
    References
    ----------
        - Tody, D. 1986, Proc. SPIE, 627, 733.
        - Tody, D. 1993, ASP Conf. Ser., 52, 173.
        - Virtanen, P., Gommers, R., Oliphant, T. E. et al. 2020,
          Nature Methods, 17, 261.
    """

    # ####################  MODULES & PACKAGES  ####################
    # ##### PANORAMA #####
    from .support.spec.spec_support import _check_fit_parameters, _check_wave_lim, _ContinuumFit, _InteractiveNormalization, _load_spectrum, _prepare_spectrum, _spectrum_file_format
    from .support.tex.tex_support import _check_output_file


    # ####################  INPUTS VALIDATION  ####################
    # ##########  VALUES  ##########
    wave, flux = _load_spectrum(spec)  # (wave, flux) pair or spectrum file
    #
    output_file = _check_output_file(output_file, allow_none=True)  # None: nothing is written
    if output_file is not None:
        _spectrum_file_format(output_file)  # Unsupported extensions are rejected before fitting
    #
    _check_wave_lim(wave_lim)
    #
    params = _check_fit_parameters(function, order, sample, naverage, low_reject, high_reject, niterate, grow, output_type)  # From 'function' to 'output_type', in this order
    #
    if not isinstance(replace, bool):
        raise TypeError(f"Replace flag 'replace' must be a boolean. Got {type(replace).__name__}.")
    #
    if not isinstance(markrej, bool):
        raise TypeError(f"Mark rejected points flag 'markrej' must be a boolean. Got {type(markrej).__name__}.")
    #
    if not isinstance(interactive, bool):
        raise TypeError(f"Interactive mode flag 'interactive' must be a boolean. Got {type(interactive).__name__}.")


    # ##########  SPECTRUM  ##########
    wave, flux = _prepare_spectrum(wave, flux, wave_lim, keep_flux_nonfin=True)  # Non-finite fluxes are kept but never used as sample points


    # ####################  COMPUTATION  ####################
    # ##########  CONTINUUM FIT  ##########
    # Spectrum, fixed parameters and resettable parameters; the sample regions are parsed against the spectrum here
    fit = _ContinuumFit(
                        wave,
                        flux,
                        function=params['function'],
                        order=params['order'],
                        sample=sample,  # Raw value: parsed by _ContinuumFit with the wavelength limits of the spectrum
                        naverage=params['naverage'],
                        low_reject=params['low_reject'],
                        high_reject=params['high_reject'],
                        niterate=params['niterate'],
                        grow=params['grow'],
                        output_type=params['output_type'],
                        replace=replace
                       )


    # ##########  OPTION 1: NON-INTERACTIVE MODE  ##########
    if not interactive:
        fit._fit()  # Automatic fit with sigma rejection; errors (e.g. too few points) are raised
        if output_file is not None:
            fit._save(output_file)  # Written with 'writespec', with a header describing the fit


    # ##########  OPTION 2: INTERACTIVE MODE  ##########
    else:
        _InteractiveNormalization(fit, output_file, markrej)._run()  # Edits 'fit' in place; returns when the window is closed


    # ####################  OUTPUT  ####################
    return wave, fit.result  # 'result' is None if no continuum could be fitted
