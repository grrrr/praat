/* Sound.cpp
 *
 * Copyright (C) 1992-2012,2014,2015,2016,2017 Paul Boersma
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This code is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this work. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * a selection of changes:
 * pb 2006/12/31 stereo
 * pb 2010/03/26 Sounds_convolve, Sounds_crossCorrelate, Sound_autocorrelate
 */

#include "Sound.h"
#include "Sound_extensions.h"
#include "NUM2.h"
#include "tensor.h"

#include "enums_getText.h"
#include "Sound_enums.h"
#include "enums_getValue.h"
#include "Sound_enums.h"

Thing_implement (Sound, Vector, 2);

autoSound Sound_clipboard;

void structSound :: v_info () {
	structDaata :: v_info ();
	const double rho_c = 400;   /* rho = 1.14 kg m-3; c = 353 m s-1; [rho c] = kg m-2 s-1 */
	double minimum = z [1] [1], maximum = minimum;
	MelderInfo_writeLine (U"Number of channels: ", ny, ny == 1 ? U" (mono)" : ny == 2 ? U" (stereo)" : U"");
	MelderInfo_writeLine (U"Time domain:");
	MelderInfo_writeLine (U"   Start time: ", xmin, U" seconds");
	MelderInfo_writeLine (U"   End time: ", xmax, U" seconds");
	MelderInfo_writeLine (U"   Total duration: ", xmax - xmin, U" seconds");
	MelderInfo_writeLine (U"Time sampling:");
	MelderInfo_writeLine (U"   Number of samples: ", nx);
	MelderInfo_writeLine (U"   Sampling period: ", dx, U" seconds");
	MelderInfo_writeLine (U"   Sampling frequency: ", Melder_single (1.0 / dx), U" Hz");
	MelderInfo_writeLine (U"   First sample centred at: ", x1, U" seconds");
	{// scope
		real80 sum = 0.0, sumOfSquares = 0.0;
		for (integer channel = 1; channel <= ny; channel ++) {
			double *amplitude = z [channel];
			for (integer i = 1; i <= nx; i ++) {
				double value = amplitude [i];
				sum += value;
				sumOfSquares += value * value;
				if (value < minimum) minimum = value;
				if (value > maximum) maximum = value;
			}
		}
		MelderInfo_writeLine (U"Amplitude:");
		MelderInfo_writeLine (U"   Minimum: ", Melder_single (minimum), U" Pascal");
		MelderInfo_writeLine (U"   Maximum: ", Melder_single (maximum), U" Pascal");
		double mean = (real) sum / (nx * ny);
		MelderInfo_writeLine (U"   Mean: ", Melder_single (mean), U" Pascal");
		MelderInfo_writeLine (U"   Root-mean-square: ", Melder_single (sqrt ((real) sumOfSquares / (nx * ny))), U" Pascal");
		double penergy = (real) sumOfSquares * dx / ny;   /* Pa2 s = kg2 m-2 s-3 */
		MelderInfo_write (U"Total energy: ", Melder_single (penergy), U" Pascal\u00B2 sec");
		double energy = penergy / rho_c;   /* kg s-2 = Joule m-2 */
		MelderInfo_writeLine (U" (energy in air: ", Melder_single (energy), U" Joule/m\u00B2)");
		double power = energy / (dx * nx);   /* kg s-3 = Watt/m2 */
		MelderInfo_write (U"Mean power (intensity) in air: ", Melder_single (power), U" Watt/m\u00B2");
		if (power != 0.0) {
			MelderInfo_writeLine (U" = ", Melder_half (10 * log10 (power / 1e-12)), U" dB");
		} else {
			MelderInfo_writeLine (U"");
		}
	}
	if (nx > 1) {
		for (integer channel = 1; channel <= ny; channel ++) {
			double stdev = stdev_scalar ({ z [channel], our nx });
			MelderInfo_writeLine (U"Standard deviation in channel ", channel, U": ", Melder_single (stdev), U" Pascal");
		}
	}
}

double structSound :: v_getMatrix (integer irow, integer icol) {
	if (irow < 1 || irow > ny) {
		if (irow == 0) {
			if (icol < 1 || icol > nx) return 0.0;
			if (ny == 1) return z [1] [icol];   // optimization
			if (ny == 2) return 0.5 * (z [1] [icol] + z [2] [icol]);   // optimization
			real80 sum = 0.0;
			for (integer channel = 1; channel <= ny; channel ++) {
				sum += z [channel] [icol];
			}
			return (real) sum / ny;
		}
		return 0.0;
	}
	if (icol < 1 || icol > nx) return 0.0;
	return z [irow] [icol];
}

double structSound :: v_getFunction2 (double x, double y) {
	integer channel = Melder_ifloor (y);
	if (channel < 0 || channel > ny || y != (double) channel) return 0.0;
	return v_getFunction1 (channel, x);
}

autoSound Sound_create (integer numberOfChannels, double xmin, double xmax, integer nx, double dx, double x1) {
	try {
		autoSound me = Thing_new (Sound);
		Matrix_init (me.get(), xmin, xmax, nx, dx, x1, 1, numberOfChannels, numberOfChannels, 1, 1);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Sound not created.");
	}
}

autoSound Sound_createSimple (integer numberOfChannels, double duration, double samplingFrequency) {
	Melder_assert (duration >= 0.0);
	Melder_assert (samplingFrequency > 0.0);
	double numberOfSamples_f = round (duration * samplingFrequency);
	if (numberOfSamples_f > (double) INT32_MAX)
		Melder_throw (U"Cannot create sounds with more than ", Melder_bigInteger (INT32_MAX), U" samples, because they cannot be saved to disk.");
	return Sound_create (numberOfChannels, 0.0, duration, (integer) (int32_t) numberOfSamples_f,
		1.0 / samplingFrequency, 0.5 / samplingFrequency);
}

/* End of file Sound.cpp */
