//   This file is part of futilities
//
//   Copyright (C) 2026 C. Ringeval
//   
//   futilities is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   futilities is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with futilities.  If not, see <https://www.gnu.org/licenses/>.


#include <stdio.h>
#include<flint.h>
#include <acb.h>
#include<acb_modular.h>

//acb_ptr ptracb;

void allocate_acb_t(int nsize, void *ptr)
{
   
  *(acb_ptr*)ptr = _acb_vec_init((slong)nsize);
  
}


int allocated_bytes_acb_t(int nsize, void *ptr)
{
  
  return _acb_vec_allocated_bytes(*(acb_srcptr*)ptr,(slong)nsize);

}


void free_acb_t(int nsize, void *ptr)
{
  _acb_vec_clear(*(acb_ptr *)ptr,(slong)nsize);
    
  return;
}


void initialize_acb_t_real(int nsize, double *x, void *ptr)
{
   
  for (int i = 0; i < nsize; i++){
    // acb_set_d(&ptracb[i],x[i]);
    //    buf = * (acb_ptr *)ptr;

    acb_set_d(*(acb_ptr *)ptr+i, x[i]);

    //printf("i= %i ",i);
    //printf("x= %.3f \n",x[i]);    
    //acb_print(*(acb_ptr *)ptr+i);
    //printf("\n");
    }
  return;
}

void initialize_acb_t_cmpx(int nsize, double *x, double *y, void *ptr)
{
  
  for (int i = 0; i < nsize; i++){
    //acb_set_d_d(&ptracb[i],x[i],y[i]);
    acb_set_d_d(*(acb_ptr *)ptr+i, x[i], y[i]);
    }
  return;
}


void get_elliptic_thetas(double theta1[2],double theta2[2], double theta3[2], double theta4[2], double z[2], double tau[2], int prec)
{

  acb_t Z,T;
  acb_t Theta1, Theta2, Theta3, Theta4;
    
  acb_init(Z);
  acb_init(T);
  
  acb_set_d_d(Z,z[0],z[1]);
  acb_set_d_d(T,tau[0],tau[1]);

  //printf("z1= %.6f \n",z[0]);
  //printf("z2= %.6f \n",z[1]);
  //printf("tau1= %.6f \n",tau[0]);
  //printf("tau2= %.6f \n",tau[1]);
  //printf("prec= %i \n",prec);
  //acb_print(Z);
  //printf("\n");
  //acb_print(T);
  //printf("\n");
  
  acb_init(Theta1);
  acb_init(Theta2);
  acb_init(Theta3);
  acb_init(Theta4);
  
  acb_modular_theta(Theta1,Theta2,Theta3,Theta4,Z,T,(slong)prec);
  
  theta1[0]= arf_get_d(arb_midref(acb_realref(Theta1)), ARF_RND_NEAR);
  theta1[1]= arf_get_d(arb_midref(acb_imagref(Theta1)), ARF_RND_NEAR);
  theta2[0]= arf_get_d(arb_midref(acb_realref(Theta2)), ARF_RND_NEAR);
  theta2[1]= arf_get_d(arb_midref(acb_imagref(Theta2)), ARF_RND_NEAR);
  theta3[0]= arf_get_d(arb_midref(acb_realref(Theta3)), ARF_RND_NEAR);
  theta3[1]= arf_get_d(arb_midref(acb_imagref(Theta3)), ARF_RND_NEAR);
  theta4[0]= arf_get_d(arb_midref(acb_realref(Theta4)), ARF_RND_NEAR);
  theta4[1]= arf_get_d(arb_midref(acb_imagref(Theta4)), ARF_RND_NEAR);
  

  acb_clear(Z);
  acb_clear(T);
  acb_clear(Theta1);
  acb_clear(Theta2);
  acb_clear(Theta3);
  acb_clear(Theta4);
  
  return;
}
