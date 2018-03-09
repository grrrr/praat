/* Matrix.cpp
 *
 * Copyright (C) 1992-2012,2013,2014,2015,2016,2017 Paul Boersma
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

#include "Matrix.h"
#include "NUM2.h"
#include "Formula.h"
#include "Eigen.h"

#include "oo_DESTROY.h"
#include "Matrix_def.h"
#include "oo_COPY.h"
#include "Matrix_def.h"
#include "oo_EQUAL.h"
#include "Matrix_def.h"
#include "oo_CAN_WRITE_AS_ENCODING.h"
#include "Matrix_def.h"
#include "oo_WRITE_TEXT.h"
#include "Matrix_def.h"
#include "oo_WRITE_BINARY.h"
#include "Matrix_def.h"
#include "oo_READ_BINARY.h"
#include "Matrix_def.h"
#include "oo_DESCRIPTION.h"
#include "Matrix_def.h"

void Matrix_drawRows (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double minimum, double maximum)
{
#ifndef NOGRAPHICS
	if (xmax <= xmin) { xmin = my xmin; xmax = my xmax; }
	if (ymax <= ymin) { ymin = my ymin; ymax = my ymax; }
	integer ixmin, ixmax, iymin, iymax;
	(void) Matrix_getWindowSamplesX (me, xmin, xmax, & ixmin, & ixmax);
	(void) Matrix_getWindowSamplesY (me, ymin, ymax, & iymin, & iymax);
	if (maximum <= minimum)
		(void) Matrix_getWindowExtrema (me, ixmin, ixmax, iymin, iymax, & minimum, & maximum);
	if (maximum <= minimum) { minimum -= 1.0; maximum += 1.0; }
	if (xmin >= xmax) return;
	Graphics_setInner (g);
	for (integer iy = iymin; iy <= iymax; iy ++) {
		Graphics_setWindow (g, xmin, xmax,
			minimum - (iy - iymin) * (maximum - minimum),
			maximum + (iymax - iy) * (maximum - minimum));
		Graphics_function (g, my z [iy], ixmin, ixmax,
			Matrix_columnToX (me, ixmin), Matrix_columnToX (me, ixmax));
	}
	Graphics_unsetInner (g);
	if (iymin < iymax)
		Graphics_setWindow (g, xmin, xmax, my y1 + (iymin - 1.5) * my dy, my y1 + (iymax - 0.5) * my dy);
#endif
}

void Matrix_drawOneContour (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double height)
{
#ifndef NOGRAPHICS
	bool xreversed = xmin > xmax, yreversed = ymin > ymax;
	if (xmax == xmin) { xmin = my xmin; xmax = my xmax; }
	if (ymax == ymin) { ymin = my ymin; ymax = my ymax; }
	if (xreversed) { double temp = xmin; xmin = xmax; xmax = temp; }
	if (yreversed) { double temp = ymin; ymin = ymax; ymax = temp; }
	integer ixmin, ixmax, iymin, iymax;
	(void) Matrix_getWindowSamplesX (me, xmin, xmax, & ixmin, & ixmax);
	(void) Matrix_getWindowSamplesY (me, ymin, ymax, & iymin, & iymax);
	if (xmin == xmax || ymin == ymax) return;
	Graphics_setInner (g);
	Graphics_setWindow (g, xreversed ? xmax : xmin, xreversed ? xmin : xmax, yreversed ? ymax : ymin, yreversed ? ymin : ymax);
	Graphics_contour (g, my z,
		ixmin, ixmax, Matrix_columnToX (me, ixmin), Matrix_columnToX (me, ixmax),
		iymin, iymax, Matrix_rowToY (me, iymin), Matrix_rowToY (me, iymax),
		height);
	Graphics_rectangle (g, xmin, xmax, ymin, ymax);
	Graphics_unsetInner (g);
#endif
}

void Matrix_drawContours (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double minimum, double maximum)
{
#ifndef NOGRAPHICS
	double border [1 + 8];
	if (xmax == xmin) { xmin = my xmin; xmax = my xmax; }
	if (ymax == ymin) { ymin = my ymin; ymax = my ymax; }
	integer ixmin, ixmax, iymin, iymax;
	(void) Matrix_getWindowSamplesX (me, xmin, xmax, & ixmin, & ixmax);
	(void) Matrix_getWindowSamplesY (me, ymin, ymax, & iymin, & iymax);
	if (maximum <= minimum)
		(void) Matrix_getWindowExtrema (me, ixmin, ixmax, iymin, iymax, & minimum, & maximum);
	if (maximum <= minimum) { minimum -= 1.0; maximum += 1.0; }
	for (integer iborder = 1; iborder <= 8; iborder ++)
		border [iborder] = minimum + iborder * (maximum - minimum) / (8 + 1);
	if (xmin == xmax || ymin == ymax) return;
	Graphics_setInner (g);
	Graphics_setWindow (g, xmin, xmax, ymin, ymax);
	Graphics_altitude (g, my z,
		ixmin, ixmax, Matrix_columnToX (me, ixmin), Matrix_columnToX (me, ixmax),
		iymin, iymax, Matrix_rowToY (me, iymin), Matrix_rowToY (me, iymax),
		8, border);
	Graphics_rectangle (g, xmin, xmax, ymin, ymax);
	Graphics_unsetInner (g);
#endif
}

void Matrix_paintContours (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double minimum, double maximum)
{
#ifndef NOGRAPHICS
	double border [1 + 30];
	if (xmax <= xmin) { xmin = my xmin; xmax = my xmax; }
	if (ymax <= ymin) { ymin = my ymin; ymax = my ymax; }
	integer ixmin, ixmax, iymin, iymax;
	(void) Matrix_getWindowSamplesX (me, xmin, xmax, & ixmin, & ixmax);
	(void) Matrix_getWindowSamplesY (me, ymin, ymax, & iymin, & iymax);
	if (maximum <= minimum)
		(void) Matrix_getWindowExtrema (me, ixmin, ixmax, iymin, iymax, & minimum, & maximum);
	if (maximum <= minimum) { minimum -= 1.0; maximum += 1.0; }
	for (integer iborder = 1; iborder <= 30; iborder ++)
		border [iborder] = minimum + iborder * (maximum - minimum) / (30 + 1);
	if (xmin >= xmax || ymin >= ymax) return;
	Graphics_setInner (g);
	Graphics_setWindow (g, xmin, xmax, ymin, ymax);
	Graphics_grey (g, my z,
		ixmin, ixmax, Matrix_columnToX (me, ixmin), Matrix_columnToX (me, ixmax),
		iymin, iymax, Matrix_rowToY (me, iymin), Matrix_rowToY (me, iymax),
		30, border);
	Graphics_rectangle (g, xmin, xmax, ymin, ymax);
	Graphics_unsetInner (g);
#endif
}

static void cellArrayOrImage (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double minimum, double maximum, bool interpolate)
{
#ifndef NOGRAPHICS
	if (xmax <= xmin) { xmin = my xmin; xmax = my xmax; }
	if (ymax <= ymin) { ymin = my ymin; ymax = my ymax; }
	integer ixmin, ixmax, iymin, iymax;
	(void) Matrix_getWindowSamplesX (me, xmin - 0.49999 * my dx, xmax + 0.49999 * my dx,
		& ixmin, & ixmax);
	(void) Matrix_getWindowSamplesY (me, ymin - 0.49999 * my dy, ymax + 0.49999 * my dy,
		& iymin, & iymax);
	if (maximum <= minimum)
		(void) Matrix_getWindowExtrema (me, ixmin, ixmax, iymin, iymax, & minimum, & maximum);
	if (maximum <= minimum) { minimum -= 1.0; maximum += 1.0; }
	if (xmin >= xmax || ymin >= ymax) return;
	Graphics_setInner (g);
	Graphics_setWindow (g, xmin, xmax, ymin, ymax);
	if (interpolate)
		Graphics_image (g, my z,
			ixmin, ixmax, Sampled_indexToX   (me, ixmin - 0.5), Sampled_indexToX   (me, ixmax + 0.5),
			iymin, iymax, SampledXY_indexToY (me, iymin - 0.5), SampledXY_indexToY (me, iymax + 0.5),
			minimum, maximum);
	else
		Graphics_cellArray (g, my z,
			ixmin, ixmax, Sampled_indexToX   (me, ixmin - 0.5), Sampled_indexToX   (me, ixmax + 0.5),
			iymin, iymax, SampledXY_indexToY (me, iymin - 0.5), SampledXY_indexToY (me, iymax + 0.5),
			minimum, maximum);
	Graphics_rectangle (g, xmin, xmax, ymin, ymax);
	Graphics_unsetInner (g);
#endif
}

void Matrix_paintImage (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double minimum, double maximum)
{
	cellArrayOrImage (me, g, xmin, xmax, ymin, ymax, minimum, maximum, true);
}

void Matrix_paintCells (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double minimum, double maximum)
{
	cellArrayOrImage (me, g, xmin, xmax, ymin, ymax, minimum, maximum, false);
}

void Matrix_paintSurface (Matrix me, Graphics g, double xmin, double xmax, double ymin, double ymax,
	double minimum, double maximum, double elevation, double azimuth)
{
#ifndef NOGRAPHICS
	if (xmax <= xmin) { xmin = my xmin; xmax = my xmax; }
	if (ymax <= ymin) { ymin = my ymin; ymax = my ymax; }
	integer ixmin, ixmax, iymin, iymax;
	(void) Matrix_getWindowSamplesX (me, xmin, xmax, & ixmin, & ixmax);
	(void) Matrix_getWindowSamplesY (me, ymin, ymax, & iymin, & iymax);
	if (maximum <= minimum)
		(void) Matrix_getWindowExtrema (me, ixmin, ixmax, iymin, iymax, & minimum, & maximum);
	if (maximum <= minimum) { minimum -= 1.0; maximum += 1.0; }
	Graphics_setInner (g);
	Graphics_setWindow (g, -1.0, 1.0, minimum, maximum);
	Graphics_surface (g, my z,
		ixmin, ixmax, Matrix_columnToX (me, ixmin), Matrix_columnToX (me, ixmax),
		iymin, iymax, Matrix_rowToY (me, iymin), Matrix_rowToY (me, iymax),
		minimum, maximum, elevation, azimuth);
	Graphics_unsetInner (g);
#endif
}

void Matrix_movie (Matrix me, Graphics g) {
#ifndef NOGRAPHICS
	autoNUMvector <double> column (1, my ny);
	double minimum = 0.0, maximum = 1.0;
	Matrix_getWindowExtrema (me, 1, my nx, 1, my ny, & minimum, & maximum);
	for (integer icol = 1; icol <= my nx; icol ++) {
		for (integer irow = 1; irow <= my ny; irow ++) {
			column [irow] = my z [irow] [icol];
		}
		Graphics_beginMovieFrame (g, & Graphics_WHITE);
		Graphics_setWindow (g, my ymin, my ymax, minimum, maximum);
		Graphics_function (g, column.peek(), 1, my ny, my ymin, my ymax);
		Graphics_endMovieFrame (g, 0.03);
	}
#endif
}

autoMatrix Matrix_readAP (MelderFile file) {
	try {
		autofile f = Melder_fopen (file, "rb");
		int16_t header [256];
		for (integer i = 0; i < 256; i ++)
			header [i] = bingeti16LE (f);
		double samplingFrequency = header [100];   // converting up (from 16 to 54 bytes)
		Melder_casual (U"Sampling frequency ", samplingFrequency);
		autoMatrix me = Matrix_create (0.0, (double) header [34], header [34] /* Number of frames. */, 1.0, 0.5,
			0.0, (double) header [35], header [35] /* Number of words per frame. */, 1.0, 0.5);
			/*Mat := MATRIX_create (Buffer.I2 [36], (* Number of words per frame. *)
							   Buffer.I2 [35], (* Number of frames. *)
							   1.0,
							   Buffer.I2 [111] / (* Samples per frame. *)
							   Buffer.I2 [101]); (* Sampling frequency. *)*/
		Melder_casual (U"... Loading ", header [34], U" frames",
			U" of ", header [35], U" words ...");
		for (integer i = 1; i <= my nx; i ++)
			for (integer j = 1; j <= my ny; j ++)
				my z [j] [i] = bingeti16LE (f);   // converting up (from 16 to 54 bytes)

		/*
		 * Get pitch frequencies.
		 */
		for (integer i = 1; i <= my nx; i ++)
			if (my z [1] [i] != 0.0)
				my z [1] [i] = - samplingFrequency / my z [1] [i];

		f.close (file);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Matrix object not read from AP file ", file);
	}
}

autoMatrix Matrix_appendRows (Matrix me, Matrix thee, ClassInfo klas) {
	try {
		autoMatrix him = Thing_newFromClass (klas).static_cast_move<structMatrix>();
		Matrix_init (him.get(), my xmin < thy xmin ? my xmin : thy xmin,
			my xmax > thy xmax ? my xmax : thy xmax,
			my nx > thy nx ? my nx : thy nx, my dx, my x1 < thy x1 ? my x1 : thy x1,
			my ymin, my ymax + (thy ymax - thy ymin), my ny + thy ny, my dy, my y1);
		for (integer irow = 1; irow <= my ny; irow ++)
			for (integer icol = 1; icol <= my nx; icol ++)
				his z [irow] [icol] = my z [irow] [icol];
		for (integer irow = 1; irow <= thy ny; irow ++)
			for (integer icol = 1; icol <= thy nx; icol ++)
				his z [irow + my ny] [icol] = thy z [irow] [icol];
		return him;
	} catch (MelderError) {
		Melder_throw (me, U" & ", thee, U": rows not appended.");
	}
}

autoMatrix Matrix_readFromRawTextFile (MelderFile file) {   // BUG: not Unicode-compatible
	try {
		autofile f = Melder_fopen (file, "rb");

		/*
		 * Count number of columns.
		 */
		integer ncol = 0;
		for (;;) {
			int kar = fgetc (f);
			if (kar == '\n' || kar == '\r' || kar == EOF) break;
			if (kar == ' ' || kar == '\t') continue;
			ncol ++;
			do {
				kar = fgetc (f);
			} while (kar != ' ' && kar != '\t' && kar != '\n' && kar != '\r' && kar != EOF);
			if (kar == '\n' || kar == '\r' || kar == EOF) break;
		}
		if (ncol == 0)
			Melder_throw (U"File empty");

		/*
		 * Count number of elements.
		 */
		rewind (f);
		integer nelements = 0;
		for (;;) {
			double element;
			if (fscanf (f, "%lf", & element) < 1) break;   // zero or end-of-file
			nelements ++;
		}

		/*
		 * Check if all columns are complete.
		 */
		if (nelements == 0 || nelements % ncol != 0)
			Melder_throw (U"The number of elements (", nelements, U") is not a multiple of the number of columns (", ncol, U").");

		/*
		 * Create simple matrix.
		 */
		integer nrow = nelements / ncol;
		autoMatrix me = Matrix_createSimple (nrow, ncol);

		/*
		 * Read elements.
		 */
		rewind (f);
		for (integer irow = 1; irow <= nrow; irow ++)
			for (integer icol = 1; icol <= ncol; icol ++)
				fscanf (f, "%lf", & my z [irow] [icol]);

		f.close (file);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Matrix object not read from raw text file ", file, U".");
	}
}

void Matrix_eigen (Matrix me, autoMatrix *out_eigenvectors, autoMatrix *out_eigenvalues) {
	try {
		if (my nx != my ny)
			Melder_throw (U"(Matrix not square.");

		autoEigen eigen = Thing_new (Eigen);
		Eigen_initFromSymmetricMatrix (eigen.get(), my z, my nx);
		autoMatrix eigenvectors = Data_copy (me);
		autoMatrix eigenvalues = Matrix_create (1.0, 1.0, 1, 1.0, 1.0, my ymin, my ymax, my ny, my dy, my y1);
		for (integer i = 1; i <= my nx; i ++) {
			eigenvalues -> z [i] [1] = eigen -> eigenvalues [i];
			for (integer j = 1; j <= my nx; j ++)
				eigenvectors -> z [i] [j] = eigen -> eigenvectors [j] [i];
		}
		*out_eigenvectors = eigenvectors.move();
		*out_eigenvalues = eigenvalues.move();
	} catch (MelderError) {
		Melder_throw (me, U": eigenstructure not computed.");
	}
}

autoMatrix Matrix_power (Matrix me, integer power) {
	try {
		if (my nx != my ny)
			Melder_throw (U"Matrix not square.");
		autoMatrix thee = Data_copy (me);
		autoMatrix him = Data_copy (me);
		for (integer ipow = 2; ipow <= power; ipow ++) {
			double **tmp = his z; his z = thy z; thy z = tmp;
			for (integer irow = 1; irow <= my ny; irow ++) {
				for (integer icol = 1; icol <= my nx; icol ++) {
					thy z [irow] [icol] = 0.0;
					for (integer i = 1; i <= my nx; i ++) {
						thy z [irow] [icol] += his z [irow] [i] * my z [i] [icol];
					}
				}
			}
		}
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": power not computed.");
	}
}

void Matrix_writeToMatrixTextFile (Matrix me, MelderFile file) {
	try {
		autofile f = Melder_fopen (file, "w");
		fprintf (f, "\"ooTextFile\"\n\"Matrix\"\n%s %s %s %s %s\n%s %s %s %s %s\n",
			Melder8_double (my xmin), Melder8_double (my xmax), Melder8_integer (my nx),
				Melder8_double (my dx), Melder8_double (my x1),
			Melder8_double (my ymin), Melder8_double (my ymax), Melder8_integer (my ny),
				Melder8_double (my dy), Melder8_double (my y1));
		for (integer i = 1; i <= my ny; i ++) {
			for (integer j = 1; j <= my nx; j ++) {
				if (j > 1) fprintf (f, " ");
				fprintf (f, "%s", Melder8_double (my z [i] [j]));
			}
			fprintf (f, "\n");
		}
		f.close (file);
	} catch (MelderError) {
		Melder_throw (me, U": not written to Matrix text file.");
	}
}

void Matrix_writeToHeaderlessSpreadsheetFile (Matrix me, MelderFile file) {
	try {
		autofile f = Melder_fopen (file, "w");
		for (integer i = 1; i <= my ny; i ++) {
			for (integer j = 1; j <= my nx; j ++) {
				if (j > 1) fprintf (f, "\t");
				fprintf (f, "%s", Melder8_single (my z [i] [j]));
			}
			fprintf (f, "\n");
		}
		f.close (file);
	} catch (MelderError) {
		Melder_throw (me, U": not saved as tab-separated file ", file);
	}
}

void Matrix_formula (Matrix me, const char32 *expression, Interpreter interpreter, Matrix target) {
	try {
		Formula_Result result;
		Formula_compile (interpreter, me, expression, kFormula_EXPRESSION_TYPE_NUMERIC, true);
		if (! target) target = me;
		for (integer irow = 1; irow <= my ny; irow ++) {
			for (integer icol = 1; icol <= my nx; icol ++) {
				Formula_run (irow, icol, & result);
				target -> z [irow] [icol] = result. numericResult;
			}
		}
	} catch (MelderError) {
		Melder_throw (me, U": formula not completed.");
	}
}

void Matrix_formula_part (Matrix me, double xmin, double xmax, double ymin, double ymax,
	const char32 *expression, Interpreter interpreter, Matrix target)
{
	try {
		if (xmax <= xmin) { xmin = my xmin; xmax = my xmax; }
		if (ymax <= ymin) { ymin = my ymin; ymax = my ymax; }
		integer ixmin, ixmax, iymin, iymax;
		(void) Matrix_getWindowSamplesX (me, xmin, xmax, & ixmin, & ixmax);
		(void) Matrix_getWindowSamplesY (me, ymin, ymax, & iymin, & iymax);
		Formula_Result result;
		Formula_compile (interpreter, me, expression, kFormula_EXPRESSION_TYPE_NUMERIC, true);
		if (! target) target = me;
		for (integer irow = iymin; irow <= iymax; irow ++) {
			for (integer icol = ixmin; icol <= ixmax; icol ++) {
				Formula_run (irow, icol, & result);
				target -> z [irow] [icol] = result. numericResult;
			}
		}
	} catch (MelderError) {
		Melder_throw (me, U": formula not completed.");
	}
}

void Matrix_scaleAbsoluteExtremum (Matrix me, double scale) {
	double extremum = 0.0;
	for (integer i = 1; i <= my ny; i ++) {
		for (integer j = 1; j <= my nx; j ++) {
			if (fabs (my z [i] [j]) > extremum) {
				extremum = fabs (my z [i] [j]);
			}
		}
	}
	if (extremum != 0.0) {
		double factor = scale / extremum;
		for (integer i = 1; i <= my ny; i ++) {
			for (integer j = 1; j <= my nx; j ++) {
				my z [i] [j] *= factor;
			}
		}
	}
}

autoMatrix TableOfReal_to_Matrix (TableOfReal me) {
	try {
		autoMatrix thee = Matrix_createSimple (my numberOfRows, my numberOfColumns);
		for (integer i = 1; i <= my numberOfRows; i ++)
			for (integer j = 1; j <= my numberOfColumns; j ++)
				thy z [i] [j] = my data [i] [j];
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": not converted to Matrix.");
	}
}

autoTableOfReal Matrix_to_TableOfReal (Matrix me) {
	try {
		autoTableOfReal thee = TableOfReal_create (my ny, my nx);
		for (integer i = 1; i <= my ny; i ++)
			for (integer j = 1; j <= my nx; j ++)
				thy data [i] [j] = my z [i] [j];
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": not converted to TableOfReal.");
	}
}

autoMatrix Table_to_Matrix (Table me) {
	try {
		autoMatrix thee = Matrix_createSimple (my rows.size, my numberOfColumns);
		for (integer icol = 1; icol <= my numberOfColumns; icol ++) {
			Table_numericize_Assert (me, icol);
		}
		for (integer irow = 1; irow <= my rows.size; irow ++) {
			TableRow row = my rows.at [irow];
			for (integer icol = 1; icol <= my numberOfColumns; icol ++) {
				thy z [irow] [icol] = row -> cells [icol]. number;
			}
		}
		return thee;
	} catch (MelderError) {
		Melder_throw (me, U": not converted to Matrix.");
	}
}

/* End of file Matrix.cpp */
