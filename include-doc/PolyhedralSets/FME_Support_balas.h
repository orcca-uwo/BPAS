#ifndef BALAS_H
#define BALAS_H

#include "FME_Support_inequality.h"
#include "FME_Support_linAlg.h"
#include "FME_Support_fme.h"

void whatIsA(inequality * , mpq_t * , int , int );

void whatIsB(inequality * , mpq_t * , int , int ,int ); 

void whatIsd(inequality * , mpq_t * , int );

void whatIsB0(mpq_t * , mpq_t * , int , int , int );

void eval(inequality , inequality * , int , int , int);

inequality * whatIsW0(mpq_t * , mpq_t * , mpq_t * , int , int , int , int );

void coefMatrix(inequality * , mpq_t * , int , int );

inequality * balasW0(inequality * , int , int , int );

inequality * projW0(inequality * , int , int ,int , int * );

int checkBoundary(mpq_t * , mpq_t * , int , int , int );






#endif
/* This file is part of the BPAS library http://www.bpaslib.org

    BPAS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BPAS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BPAS.  If not, see <http://www.gnu.org/licenses/>.

    Copyright:
        Mohammadali Asadi <masadi4@uwo.ca>
        Alexander Brandt <abrandt5@uwo.ca>
        Changbo Chen <changbo.chen@hotmail.com>
        Svyatoslav Covanov <svyatoslav.covanov@loria.fr>
        Farnam Mansouri <mansouri.farnam@gmail.com>
        Davood Mohajerani <mohajerani.d@gmail.com>
        Robert Moir <robert@moir.net>
        Marc Moreno Maza  <moreno@csd.uwo.ca>
        Delaram Talaashrafi <dtalaash@uwo.ca>
        Amha Tsegaye <atsegaye@uwo.ca>
        Linxiao Wang <lwang739@uwo.ca>
        Ning Xie <nxie6@csd.uwo.ca>
        Yuzhen Xie <yuzhenxie@yahoo.ca>

*/


