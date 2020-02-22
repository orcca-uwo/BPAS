#ifndef _SMALLPRIMEFIELD_SUPPORT_H_
#define _SMALLPRIMEFIELD_SUPPORT_H_

#ifdef __cplusplus
    extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

long long int smallprimefield_covert_in(long long int _a, long int prime, long long int R);
long long int smallprimefield_covert_out(long long int a, long int prime, long long int R, long long int Pp);
long long int smallprimefield_getPp(long int prime, long long int R);
long long int smallprimefield_add(long long int a, long long int b, long int prime);
long long int smallprimefield_sub( long long int a, long long int b, long int prime);
long long int smallprimefield_multi(long long int a, long long int b, long int prime, long long int R, long long int Pp);
long long int smallprimefield_inverse(long long int a, long long int prime, long long int R, long long int Pp);

#ifdef __cplusplus
}
#endif

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


