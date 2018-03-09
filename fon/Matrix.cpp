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

Thing_implement (Matrix, SampledXY, 2);

void structMatrix :: v_info () {
	structDaata :: v_info ();
	double minimum = 0.0, maximum = 0.0;
	Matrix_getWindowExtrema (this, 1, our nx, 1, our ny, & minimum, & maximum);
	MelderInfo_writeLine (U"xmin: ", our xmin);
	MelderInfo_writeLine (U"xmax: ", our xmax);
	MelderInfo_writeLine (U"Number of columns: ", our nx);
	MelderInfo_writeLine (U"dx: ", our dx, U" (-> sampling rate ", 1.0 / our dx, U" )");
	MelderInfo_writeLine (U"x1: ", our x1);
	MelderInfo_writeLine (U"ymin: ", our ymin);
	MelderInfo_writeLine (U"ymax: ", our ymax);
	MelderInfo_writeLine (U"Number of rows: ", our ny);
	MelderInfo_writeLine (U"dy: ", our dy, U" (-> sampling rate ", 1.0 / our dy, U" )");
	MelderInfo_writeLine (U"y1: ", our y1);
	MelderInfo_writeLine (U"Minimum value: ", minimum);
	MelderInfo_writeLine (U"Maximum value: ", maximum);
}

void structMatrix :: v_readText (MelderReadText text, int formatVersion) {
	if (formatVersion < 0) {
		our xmin = texgetr64 (text);
		our xmax = texgetr64 (text);
		our ymin = texgetr64 (text);
		our ymax = texgetr64 (text);
		our nx = texgeti32 (text);
		our ny = texgeti32 (text);
		our dx = texgetr64 (text);
		our dy = texgetr64 (text);
		our x1 = texgetr64 (text);
		our y1 = texgetr64 (text);
	} else {
		Matrix_Parent :: v_readText (text, formatVersion);
	}
	if (our xmin > our xmax)
		Melder_throw (U"xmin should be less than or equal to xmax.");
	if (our ymin > our ymax)
		Melder_throw (U"ymin should be less than or equal to ymax.");
	if (our nx < 1)
		Melder_throw (U"nx should be at least 1.");
	if (our ny < 1)
		Melder_throw (U"ny should be at least 1.");
	if (our dx <= 0.0)
		Melder_throw (U"dx should be greater than 0.0.");
	if (our dy <= 0.0)
		Melder_throw (U"dy should be greater than 0.0.");
	our z = NUMmatrix_readText_r64 (1, our ny, 1, our nx, text, "z");
}

double structMatrix :: v_getValueAtSample (integer isamp, integer ilevel, int unit) {
	double value = our z [ilevel] [isamp];
	return ( isdefined (value) ? our v_convertStandardToSpecialUnit (value, ilevel, unit) : undefined );
}

double structMatrix :: v_getMatrix (integer irow, integer icol) {
	if (irow < 1 || irow > our ny) return 0.0;
	if (icol < 1 || icol > our nx) return 0.0;
	return z [irow] [icol];
}

double structMatrix :: v_getFunction2 (double x, double y) {
	double rrow = (y - our y1) / our dy + 1.0;
	double rcol = (x - our x1) / our dx + 1.0;
	integer irow = Melder_ifloor (rrow), icol = Melder_ifloor (rcol);
	double drow = rrow - irow, dcol = rcol - icol;
	double z1 = irow < 1 || irow >  our ny || icol < 1 || icol >  our nx ? 0.0 : z [irow]     [icol];
	double z2 = irow < 0 || irow >= our ny || icol < 1 || icol >  our nx ? 0.0 : z [irow + 1] [icol];
	double z3 = irow < 1 || irow >  our ny || icol < 0 || icol >= our nx ? 0.0 : z [irow]     [icol + 1];
	double z4 = irow < 0 || irow >= our ny || icol < 0 || icol >= our nx ? 0.0 : z [irow + 1] [icol + 1];
	return (1.0 - drow) * (1.0 - dcol) * z1 + drow * (1.0 - dcol) * z2 + (1.0 - drow) * dcol * z3 + drow * dcol * z4;
}

autoMatrix Matrix_create
	(double xmin, double xmax, integer nx, double dx, double x1,
	 double ymin, double ymax, integer ny, double dy, double y1)
{
	try {
		autoMatrix me = Thing_new (Matrix);
		Matrix_init (me.get(), xmin, xmax, nx, dx, x1, ymin, ymax, ny, dy, y1);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Matrix object not created.");
	}
}

autoMatrix Matrix_createSimple (integer numberOfRows, integer numberOfColumns) {
	try {
		autoMatrix me = Thing_new (Matrix);
		Matrix_init (me.get(), 0.5, numberOfColumns + 0.5, numberOfColumns, 1, 1,
			0.5, numberOfRows + 0.5, numberOfRows, 1, 1);
		return me;
	} catch (MelderError) {
		Melder_throw (U"Matrix object not created.");
	}
}

double Matrix_columnToX (Matrix me, double column) { return my x1 + (column - 1.0) * my dx; }   // FIXME inline and use Sampled

double Matrix_rowToY (Matrix me, double row) { return my y1 + (row - 1.0) * my dy; }

double Matrix_xToColumn (Matrix me, double x) { return (x - my x1) / my dx + 1.0; }

integer Matrix_xToLowColumn (Matrix me, double x) { return Melder_ifloor (Matrix_xToColumn (me, x)); }

integer Matrix_xToHighColumn (Matrix me, double x) { return Melder_iceiling (Matrix_xToColumn (me, x)); }

integer Matrix_xToNearestColumn (Matrix me, double x) { return Melder_iround (Matrix_xToColumn (me, x)); }

double Matrix_yToRow (Matrix me, double y) { return (y - my y1) / my dy + 1.0; }

integer Matrix_yToLowRow (Matrix me, double y) { return Melder_ifloor (Matrix_yToRow (me, y)); }

integer Matrix_yToHighRow (Matrix me, double y) { return Melder_iceiling (Matrix_yToRow (me, y)); }

integer Matrix_yToNearestRow (Matrix me, double y) { return Melder_iround (Matrix_yToRow (me, y)); }

integer Matrix_getWindowSamplesX (Matrix me, double xmin, double xmax, integer *ixmin, integer *ixmax) {
	*ixmin = 1 + Melder_iceiling ((xmin - my x1) / my dx);
	*ixmax = 1 + Melder_ifloor   ((xmax - my x1) / my dx);
	if (*ixmin < 1) *ixmin = 1;
	if (*ixmax > my nx) *ixmax = my nx;
	if (*ixmin > *ixmax) return 0;
	return *ixmax - *ixmin + 1;
}

integer Matrix_getWindowSamplesY (Matrix me, double ymin, double ymax, integer *iymin, integer *iymax) {
	*iymin = 1 + Melder_iceiling ((ymin - my y1) / my dy);
	*iymax = 1 + Melder_ifloor   ((ymax - my y1) / my dy);
	if (*iymin < 1) *iymin = 1;
	if (*iymax > my ny) *iymax = my ny;
	if (*iymin > *iymax) return 0;
	return *iymax - *iymin + 1;
}

integer Matrix_getWindowExtrema (Matrix me, integer ixmin, integer ixmax, integer iymin, integer iymax,
	double *minimum, double *maximum)
{
	if (ixmin == 0) ixmin = 1;
	if (ixmax == 0) ixmax = my nx;
	if (iymin == 0) iymin = 1;
	if (iymax == 0) iymax = my ny;
	if (ixmin > ixmax || iymin > iymax) return 0;
	*minimum = *maximum = my z [iymin] [ixmin];
	for (integer iy = iymin; iy <= iymax; iy ++) {
		for (integer ix = ixmin; ix <= ixmax; ix ++) {
			if (my z [iy] [ix] < *minimum) *minimum = my z [iy] [ix];
			if (my z [iy] [ix] > *maximum) *maximum = my z [iy] [ix];
		}
	}
	return (ixmax - ixmin + 1) * (iymax - iymin + 1);
}

double Matrix_getValueAtXY (Matrix me, double x, double y) {
	real row_real = (y - my y1) / my dy + 1.0;
	real col_real = (x - my x1) / my dx + 1.0;
	/*
	 * We imagine a unit square around every (xi, yi) point in the matrix.
	 * For (x, y) values outside the union of these squares, the z value is undefined.
	 */
	if (row_real < 0.5 || row_real > my ny + 0.5) return undefined;
	if (col_real < 0.5 || col_real > my nx + 0.5) return undefined;
	/*
	 * Determine the four nearest (xi, yi) points.
	 */
	integer bottomRow = Melder_ifloor (row_real);   // 0 <= bottomRow <= my ny
	integer topRow = bottomRow + 1;         // 1 <= topRow <= my ny + 1
	integer leftCol = Melder_ifloor (col_real);     // 0 <= leftCol <= my nx
	integer rightCol = leftCol + 1;         // 1 <= rightCol <= my nx + 1
	real drow = row_real - bottomRow;    // 0.0 <= drow < 1.0
	real dcol = col_real - leftCol;      // 0.0 <= dcol < 1.0
	/*
	 * If adjacent points exist
	 * (i.e., both row numbers are between 1 and my ny,
	 *  or both column numbers are between 1 and my nx),
	 * we do linear interpolation.
	 * If not, we do constant extrapolation,
	 * which can be simulated by an interpolation between equal z values.
	 */
	if (bottomRow < 1) bottomRow = 1;         // 1 <= bottomRow <= my ny
	if (topRow > my ny) topRow = my ny;       // 1 <= topRow <= my ny
	if (leftCol < 1) leftCol = 1;             // 1 <= leftCol <= my nx
	if (rightCol > my nx) rightCol = my nx;   // 1 <= rightCol <= my nx
	return (1.0 - drow) * (1.0 - dcol) * my z [bottomRow] [leftCol] +
		drow * (1.0 - dcol) * my z [topRow] [leftCol] +
		(1.0 - drow) * dcol * my z [bottomRow] [rightCol] +
		drow * dcol * my z [topRow] [rightCol];
}

double Matrix_getSum (Matrix me) {
	real80 sum = 0.0;
	for (integer irow = 1; irow <= my ny; irow ++)
		for (integer icol = 1; icol <= my nx; icol ++)
			sum += my z [irow] [icol];
	return (real) sum;
}

double Matrix_getNorm (Matrix me) {
	real80 sum = 0.0;
	for (integer irow = 1; irow <= my ny; irow ++)
		for (integer icol = 1; icol <= my nx; icol ++)
			sum += my z [irow] [icol] * my z [irow] [icol];
	return sqrt ((real) sum);
}

/* End of file Matrix.cpp */
