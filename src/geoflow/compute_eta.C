/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: step.C 164 2014-07-10 12:11:22 haghakha $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

double compute_eta(HashTable* El_Table, PileProps *pileprops) {

	double local_eta = 0., eta = 0., *dx, phi = 0.0, height = 0.;
	HashEntryPtr* buck = El_Table->getbucketptr();

	for (int i = 0; i < El_Table->get_no_of_buckets(); i++)
		if (*(buck + i)) {
			HashEntryPtr currentPtr = *(buck + i);
			while (currentPtr) {
				Element* Curr_El = (Element*) (currentPtr->value);
				if (Curr_El->get_adapted_flag() > 0) {

					dx = Curr_El->get_dx();
					phi = *(Curr_El->get_state_vars());
					height = *(Curr_El->get_state_vars() + 1);

					local_eta += height*dx[0] * dx[1] * phi * (phi * phi - 1);
				}
				currentPtr = currentPtr->next;
			}
		}
	//cout<<"local eta is  "<<local_eta<<"  and volume is  "<<pileprops->pilevol<<endl;

	local_eta /=  pileprops->pilevol;

//	cout<<"local eta is  "<<local_eta<<"  and volume is  "<<pileprops->pilevol<<endl;
	MPI_Allreduce(&local_eta, &eta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return eta;
}
