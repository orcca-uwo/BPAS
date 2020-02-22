
/**
 * A simple timer for calculating user time.
 */

#ifndef _UNIX_TIMER_H_
#define _UNIX_TIMER_H_

#include <sys/time.h>
#include <sys/resource.h>

typedef struct timeval timer_time;
typedef struct {
	timer_time start_time;
} timer_id;

timer_time user_time() {
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);

	return usage.ru_utime;
}

timer_id start_timer() {
	timer_time start_time = user_time();
	timer_id id;
	id.start_time = start_time;
	return id;
}

timer_time elapsed_time(timer_id* timerID) {
	if (timerID == NULL) {
		return user_time();
	}

	timer_time t = user_time();
	t.tv_sec -= timerID->start_time.tv_sec;
	t.tv_usec -= timerID->start_time.tv_usec; 
	return t;
}


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


