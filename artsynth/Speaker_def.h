/* Speaker_def.h
 *
 * Copyright (C) 1992-2005,2011,2015-2017 Paul Boersma
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


#define ooSTRUCT Speaker_CordDimensions
oo_DEFINE_STRUCT (Speaker_CordDimensions)

	oo_INT16 (numberOfMasses)
	oo_FLOATTYPE (length)

oo_END_STRUCT (Speaker_CordDimensions)
#undef ooSTRUCT


#define ooSTRUCT Speaker_CordSpring
oo_DEFINE_STRUCT (Speaker_CordSpring)

	oo_FLOATTYPE (thickness)
	oo_FLOATTYPE (mass)
	oo_FLOATTYPE (k1)

oo_END_STRUCT (Speaker_CordSpring)
#undef ooSTRUCT


#define ooSTRUCT Speaker_GlottalShunt
oo_DEFINE_STRUCT (Speaker_GlottalShunt)

	oo_FLOATTYPE (Dx)
	oo_FLOATTYPE (Dy)
	oo_FLOATTYPE (Dz)

oo_END_STRUCT (Speaker_GlottalShunt)
#undef ooSTRUCT


#define ooSTRUCT Speaker_Velum
oo_DEFINE_STRUCT (Speaker_Velum)   // V

	oo_FLOATTYPE (x)
	oo_FLOATTYPE (y)
	oo_FLOATTYPE (a)

oo_END_STRUCT (Speaker_Velum)
#undef ooSTRUCT


#define ooSTRUCT Speaker_Palate
oo_DEFINE_STRUCT (Speaker_Palate)   // OM

	oo_FLOATTYPE (radius)

oo_END_STRUCT (Speaker_Palate)
#undef ooSTRUCT


#define ooSTRUCT Speaker_Tip
oo_DEFINE_STRUCT (Speaker_Tip)

	oo_FLOATTYPE (length)

oo_END_STRUCT (Speaker_Tip)
#undef ooSTRUCT


#define ooSTRUCT Speaker_Alveoli
oo_DEFINE_STRUCT (Speaker_Alveoli)

	oo_FLOATTYPE (x)
	oo_FLOATTYPE (y)
	oo_FLOATTYPE (a)

oo_END_STRUCT (Speaker_Alveoli)
#undef ooSTRUCT


#define ooSTRUCT Speaker_TeethCavity
oo_DEFINE_STRUCT (Speaker_TeethCavity)

	oo_FLOATTYPE (dx1)
	oo_FLOATTYPE (dx2)
	oo_FLOATTYPE (dy)

oo_END_STRUCT (Speaker_TeethCavity)
#undef ooSTRUCT


#define ooSTRUCT Speaker_LowerTeeth
oo_DEFINE_STRUCT (Speaker_LowerTeeth)   // rest position of J

	oo_FLOATTYPE (r)
	oo_FLOATTYPE (a)

oo_END_STRUCT (Speaker_LowerTeeth)
#undef ooSTRUCT


#define ooSTRUCT Speaker_UpperTeeth
oo_DEFINE_STRUCT (Speaker_UpperTeeth)   // U

	oo_FLOATTYPE (x)
	oo_FLOATTYPE (y)

oo_END_STRUCT (Speaker_UpperTeeth)
#undef ooSTRUCT


#define ooSTRUCT Speaker_Lip
oo_DEFINE_STRUCT (Speaker_Lip)

	oo_FLOATTYPE (dx)
	oo_FLOATTYPE (dy)

oo_END_STRUCT (Speaker_Lip)
#undef ooSTRUCT


#define ooSTRUCT Speaker_Nose
oo_DEFINE_STRUCT (Speaker_Nose)

	oo_FLOATTYPE (Dx)
	oo_FLOATTYPE (Dz)
	oo_FLOATTYPE_ARRAY (weq, 14, 14)

oo_END_STRUCT (Speaker_Nose)
#undef ooSTRUCT


#define ooSTRUCT Speaker
oo_DEFINE_CLASS (Speaker, Daata)

	oo_FLOATTYPE (relativeSize)   // different for female, male, child

	/* In the larynx. */

	oo_STRUCT (Speaker_CordDimensions, cord)
	oo_STRUCT (Speaker_CordSpring, lowerCord)
	oo_STRUCT (Speaker_CordSpring, upperCord)
	oo_STRUCT (Speaker_GlottalShunt, shunt)

	/* Above the larynx. */

	oo_STRUCT (Speaker_Velum, velum)
	oo_STRUCT (Speaker_Palate, palate)
	oo_STRUCT (Speaker_Tip, tip)
	oo_FLOATTYPE (neutralBodyDistance)
	oo_STRUCT (Speaker_Alveoli, alveoli)
	oo_STRUCT (Speaker_TeethCavity, teethCavity)
	oo_STRUCT (Speaker_LowerTeeth, lowerTeeth)
	oo_STRUCT (Speaker_UpperTeeth, upperTeeth)
	oo_STRUCT (Speaker_Lip, lowerLip)
	oo_STRUCT (Speaker_Lip, upperLip)

	/* In the nasal cavity. */

	oo_STRUCT (Speaker_Nose, nose)

oo_END_CLASS (Speaker)
#undef ooSTRUCT

/* End of file Speaker_def.h */
