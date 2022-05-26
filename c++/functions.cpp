//rego: Automatic time series forecasting and missing value imputation.
//
//Copyright (C) Davide Altomare and David Loris <https://channelattribution.io>
//
//This source code is licensed under the MIT license found in the
//LICENSE file in the root directory of this source tree. 

#define language_cpp 
//#define language_python
//#define language_R

#include <iostream>
#include <vector>
#include <set>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <sstream>
#include <list>
#include <string>
#include <random>
#include <numeric>
#include <time.h> 
#include <thread>
#include <map>
#include <algorithm>
#include <list>
#include <limits> 
#include <functional>
#include <ctime>


#define OPTIM_ENABLE_ARMA_WRAPPERS
#ifdef language_py
  #ifdef _WIN32
    #define ARMA_DONT_USE_LAPACK
    #define ARMA_DONT_USE_BLAS
  #endif
#endif
#define ARMA_USE_CXX11
#define ARMA_64BIT_WORD
#define ARMA_DONT_PRINT_ERRORS

#include <armadillo>
#include <optim.hpp>

#ifdef language_R
 #define __GXX_EXPERIMENTAL_CXX0X__ 1

 #include <Rcpp.h>
 //#include <RcppArmadillo.h>
 //#define USE_RCPP_ARMADILLO
 
 #ifndef BEGIN_RCPP
 #define BEGIN_RCPP
 #endif
  
 #ifndef END_RCPP
 #define END_RCPP
 #endif
 
 using namespace Rcpp;
#endif


#ifdef language_py
 #include <Python.h>
#endif


using namespace std;
using namespace arma;


#define uli unsigned long int
using svec1=vector<double>;
using svec2=vector< vector <double> >;
using svec3=vector < vector< vector <double> > >;
using svec4=vector < vector < vector< vector <double> > > >;


//------------------------------------------------------------------------------------------------------------------------
//GENERAL FUNCTIONS
//------------------------------------------------------------------------------------------------------------------------


double lfactorial(uli n)
{
  
 double x;

 x=lgamma(n+1);

 return(x);

} //end function


double lchoose(uli n, uli k)
{
  
  double x;
  
  x=lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1);

  return(x);

} //end function


mat sub_mat(mat* M, vec vr, vec vc)
{
  uli r, c, lvr, lvc;
  lvr=vr.n_elem;
  lvc=vc.n_elem;
  mat Q(lvr,lvc);

  for(c=0; c<lvc; c++){
   for(r=0; r<lvr; r++){ 
    Q(r,c)=(*M)(vr(r),vc(c)); 
   }
  }

  return Q;

} //end function


mat pow_vec(vec* v, vec* w)
{

   uli k, lv; 
   vec vv;

   lv=(*v).n_elem;
   vv=zeros<vec>(lv);

   for(k=0; k<lv; k++){
    vv(k)=pow((*v)(k),(*w)(k));
   }

 return vv;

}

mat Cholesky(mat* M)
{
    double n=(*M).n_rows;
    mat lower(n,n);
 
    // Decomposing a matrix into Lower Triangular
    for (uli i = 0; i < n; i++) {
        for (uli j = 0; j <= i; j++) {
            uli sum = 0;
 
            if (j == i) // summation for diagonals
            {
                for (uli k = 0; k < j; k++){
                    sum += pow(lower(j,k), 2);
                    lower(j,j) = sqrt((*M)(j,j) - sum);
                }
            } else {
 
                // Evaluating L(i, j) using L(j, j)
                for (uli k = 0; k < j; k++){
                    sum += (lower(i,k) * lower(j,k));
                    lower(i,j) = ((*M)(i,j) - sum) / lower(j,j);
                }
            }
        }
    }

    return(lower);
 
}


vec forward_sub(mat* L, vec* b)
{
    /*x = forward_sub(L, b) is the solution to L x = b
       L must be a lower-triangular matrix
       b must be a vector of the same leading dimension as L
    */
    
    double n = (*L).n_cols;
    vec x = zeros<vec>(n);
    double tmp;
    for(double i=0; i<n; ++i){
        tmp = (*b)(i);
         for(double j=0; j<i; ++j){
          tmp -= (*L)(i,j) * x(j);
         }

        x(i) = tmp / (*L)(i,i);
    }
    return(x);
}

vec back_sub(mat* U, vec* b)
{
    /*x = back_sub(U, b) is the solution to U x = b
       U must be an upper-triangular matrix
       b must be a vector of the same leading dimension as U
    */
    
    double n = (*U).n_cols; //it must be double
    vec x = zeros<vec>(n);
    double tmp;
    for(double i=n-1; i>=0; i--){
      tmp = (*b)(i);
       for(uli j=i; j<n; ++j){
        tmp -= (*U)(i,j) * x(j);
        x(i) = tmp / (*U)(i,i);
       }
    }

    return(x);

}

vec solve0(mat* A, vec* b){
 
  mat L=Cholesky(A);
  mat Lt=trans(L);

  vec vy=forward_sub(&L,b);
  vec vx=back_sub(&Lt, &vy);
  
  return(vx);

}

vec solve_linear(mat* A, vec* b){
  
 #ifdef ARMA_DONT_USE_LAPACK
  return(solve0(A,b));
 #else
  return(solve(*A,*b));
 #endif
 

}


//------------------------------------------------------------------------------------------------------------------------
//FUNCTIONS FOR STORING MODELS IN A BINARY TREE
//------------------------------------------------------------------------------------------------------------------------

//this function add a model "M" to "tree"

field<mat> add_to_tree(vec* M, double lM, uli nM, mat* tree, double ltree)
{

 uli j,k;

 field<mat> Res(2,1);

 if(nM==0){
    
  for(j=0; j<=lM; j++){
 
   if((*M)(j)==1){(*tree)(j,0)=j+1;}
   else{(*tree)(j,1)=j+1;}
 
  }
 
  ltree=lM; 
  (*tree).row(ltree)=(*tree).row(ltree)*0+nM;

 }

 uli z,h,iM=0;

 if(nM>0){ //if1
   
  z=0;
  h=ltree+1;
  
  for(j=0; j<=lM; j++){ //for1  
  
   iM=1-(*M)(j);
      
   if(!is_finite((*tree)(z,iM)) && (j<=lM) ){ //if2
     
    (*tree)(z,iM)=h;
     
    for(k=(j+1); k<=lM; k++){ //for2
  
      if((*M)(k)==1){(*tree)(h,0)=h+1;} else{(*tree)(h,1)=h+1;}
      
      h=h+1;
    
    } //end for2    

    iM=1-(*M)(lM);
    ltree=h-1;
    break;
   
   } //end if2
 
   if(j==lM){(*tree)(z,iM)=nM; ltree=ltree+1; break;}
  
   if((*tree)(z,iM)>=0){z=(*tree)(z,iM);}
  
  } //end for1
  

  (*tree)(ltree,iM)=(*tree)(ltree,iM)*0+nM;

 } //end if1


 Res(0,0)=(*tree);
 Res(1,0)=ltree;

 return Res;

} //end function


//this function returns the possible movements from a model "M", given all the previous models visited and stored in "tree"

vec mov_tree(mat* tree, vec* M, uli lM, vec* vlM, uli max_lM)
{

 uli q, k, z, h, iM2;
 double sumM;
 vec mov(lM+1); 
 uvec imov, umov; 
 vec mov2, M2;
 
 
 mov.fill(-1);
 sumM=sum(*M);
 q=0;

 for(k=0; k<=lM; k++){ //for1
  
  M2=(*M); 
  M2(k)=1-(*M)(k);
  z=0;
  
  for(h=0; h<=lM; h++){ //for2
    
   iM2=1-M2(h);
   if(!is_finite((*tree)(z,iM2))){mov(q)=k; q=q+1; break;} else{z=(*tree)(z,iM2);}
   
    } //end for2
  
 
 } //end for1

 imov=find(mov>-1);
 
 if(!imov.is_empty()){
 
  mov=mov.elem(imov);
  umov=conv_to<uvec>::from(mov);

  if(sumM>=max_lM){
  
   mov2=zeros<vec>(lM+1);
   mov2.elem(umov)=ones<vec>(mov.n_elem);
   mov=(mov2%(*M))%(*vlM);
   imov=find(mov>0);
   if(!imov.is_empty()){mov=mov.elem(imov); mov=mov-1;} else{mov=datum::nan;}
  }

 } else {mov=datum::nan;}
 

 return mov;

} //end function



//------------------------------------------------------------------------------------------------------------------------
//FUNCTIONS FOR BAYESIAN STOCHASTIC SEARCH
//------------------------------------------------------------------------------------------------------------------------


double log_H_h_i(double mu, double sigma, double h, double i)
{

 double x;

 x=lfactorial(2*h)+i*log(sigma)-lfactorial(i)+(2*h-2*i)*log(abs(mu))-lfactorial(2*h-2*i);

 return x;

} //end function



double log_FBF_Ga_Gb(vec* G_a, vec* G_b, uli edge, mat* edges, mat* YtY, uli add, double n, double h)
{
  
  uli e1, e2, iwi;
  double i, p, b, S2, mu, sigma, logS2, ilogS2, logHhi, ilog4, log_num1, log_den1, log_w_1, log_num0i0, log_den0i0, log_w_0, log_FBF_unpasso;
  vec V1, V2, G1, V11, pa1, pa0, betah, vv(1), z1;
  uvec iw, ipa1;
  mat e, yty, XtX, invXtX;   
  vec Xty;
  
  e=(*edges).row(edge);
  e1=e(0);
  e2=e(1);

  V1=(*edges).col(0);
  V2=(*edges).col(1);
    
  if(add==1){G1=(*G_a);}else{G1=(*G_b);}

  V11=(V1+1)%G1;
  iw=find(V2==e2); pa1=V11.elem(iw);
  iw=find(pa1>0); pa1=pa1.elem(iw); pa1=pa1-1;
 
  iw=find(pa1!=e1); if(!iw.is_empty()){pa0=pa1.elem(iw);}else{pa0=datum::nan;}
 
  p=pa1.n_elem;
  b=(p+2*h+1)/n;

  yty=(*YtY)(e2,e2);
    
  // //calcolo w1
   
  vv(0)=e2; Xty=conv_to<vec>::from(sub_mat(YtY,pa1,vv));
  XtX=sub_mat(YtY,pa1,pa1);
  betah=solve_linear(&XtX,&Xty);

  S2=conv_to<double>::from(yty-(trans(Xty)*betah));
  
  iw=find(pa1==e1); mu=conv_to<double>::from(betah.elem(iw));
  iwi=conv_to<uli>::from(iw); 
  z1=zeros<vec>(pa1.n_elem);
  z1(iwi)=1;
  z1=solve_linear(&XtX,&z1);
  
  sigma=conv_to<double>::from(z1.elem(iw));
  
  if(S2>0){
   log_w_1=(-n*(1-b)/2)*log(datum::pi*b*S2);
   logS2=log(S2);
   log_num1=-datum::inf;
   log_den1=-datum::inf;

   for(i=0; i<=h; i++){
   
    ilogS2=i*logS2;
    logHhi=log_H_h_i(mu,sigma,h,i);
    ilog4=-i*log(4);
      
    log_num1=log_add(log_num1, (ilog4+logHhi+lgamma((n-p-2*i)/2)+ilogS2));
    log_den1=log_add(log_den1, (ilog4+logHhi+lgamma((n*b-p-2*i)/2)+ilogS2));

   }
    
   log_w_1=log_w_1+log_num1-log_den1;
  }else{
   log_w_1=datum::inf;
  }
     
  //calcolo w0

  if(!pa0.is_finite()){p=0;}else{p=pa0.n_elem;}

  log_num0i0=lgamma((n-p)/2);
  log_den0i0=lgamma((n*b-p)/2);

  if(p==0){S2=conv_to<double>::from(yty);}
  else{
   vv(0)=e2; Xty=conv_to<vec>::from(sub_mat(YtY,pa0,vv));
   XtX=sub_mat(YtY,pa0,pa0);
   betah=solve_linear(&XtX,&Xty); 

   S2=conv_to<double>::from(yty-(trans(Xty)*betah));
  }

  if(S2>0){
   log_w_0=(-(n*(1-b)/2))*log(datum::pi*b*S2)+log_num0i0-log_den0i0;
  }else{
   log_w_0=datum::inf;
  }
  
  //calcolo FBF

  if(add==1){log_FBF_unpasso=log_w_1-log_w_0;}
  else{log_FBF_unpasso=log_w_0-log_w_1;} 
  
  if(!is_finite(log_FBF_unpasso)){
   log_FBF_unpasso=0;
  }

  return log_FBF_unpasso;

} // end function




field<mat> FBF_heart(double nt, mat* YtY, vec* vG_base, double lcv, vec* vlcv, mat* edges, double n_tot_mod, double C, double maxne, double h, bool univariate)
{
    
   uli t, add, edge, imq, limodR, s;
   double ltree, lM, sum_log_FBF, log_FBF_G, log_pi_G, log_num_MP_G, sum_log_RSMP, n_mod_r, log_FBF_t, log_FBF1;
   vec M_log_FBF, log_num_MP, log_sume, G, imod_R, M_log_RSMP, pRSMP, mov, vlM, qh, G_t, M_q, M_P;
   uvec iw;
   mat tree, SM, M_G; 
   field<mat> treeRes, Res(4,1);
   uword i_n_mod_r, imaxe;
 
   M_G=zeros<mat>(lcv,n_tot_mod);  
   M_P=zeros<vec>(n_tot_mod); 
   M_log_FBF=zeros<vec>(n_tot_mod); 
   log_num_MP=zeros<vec>(n_tot_mod); 
   M_q=zeros<vec>(lcv); 
   tree=zeros<mat>(n_tot_mod*lcv,2); tree.fill(datum::nan);
   ltree=datum::nan;
   lM=lcv-1; 
 
   sum_log_FBF=-datum::inf; 
   log_sume=zeros<vec>(lcv); log_sume.fill(-datum::inf);
  
   M_log_RSMP=zeros<vec>(n_tot_mod);
   sum_log_RSMP=-datum::inf;
   imod_R=zeros<vec>(n_tot_mod);

   lM=lcv-1;

   Col<uli> vexit(1);

   for(t=0; t<lcv; t++){ //for1
       
    G=(*vG_base);
    G(t)=1-(*vG_base)(t);
    add=G(t);
    edge=t;

    log_FBF_G=log_FBF_Ga_Gb(&G,vG_base,edge,edges,YtY,add,nt,h);
    
    M_G.col(t)=G;
    
    treeRes=add_to_tree(&G,lM,t,&tree,ltree);
    tree=treeRes(0,0);
    ltree=conv_to<double>::from(treeRes(1,0));
        
    M_log_FBF(t)=log_FBF_G;
    log_pi_G=-log(lcv+1)-lchoose(lcv,sum(G));
    log_num_MP_G=log_FBF_G+log_pi_G;
    log_num_MP(t)=log_num_MP_G;

    sum_log_FBF=log_add(sum_log_FBF, log_num_MP_G);
  
    for(imq=0; imq<lcv; imq++){
     if(G(imq)==1){log_sume(imq)=log_add(log_sume(imq), log_num_MP_G);}
    }
    
    M_q=exp(log_sume-sum_log_FBF);
   
    M_log_RSMP(t)=log_num_MP_G;
    sum_log_RSMP=log_add(sum_log_RSMP, log_num_MP_G);
   
    imod_R(t)=t;
 
   } //end for1

   vec vtmp1,vtmp2;

   if(univariate==0){
   
    limodR=t-1; 
    s=lcv;

    
    while(t<n_tot_mod){ //while1

     pRSMP=exp(M_log_RSMP.subvec(0,limodR)-sum_log_RSMP);
	   pRSMP.max(i_n_mod_r);
     
     n_mod_r=imod_R(i_n_mod_r);
         
     G=M_G.col(n_mod_r);
     G_t=G;
     log_FBF_t=M_log_FBF(n_mod_r);
     
     
     vlM=(*vlcv)+1;
     mov=mov_tree(&tree,&G,lM,&vlM,maxne);
   
     if(!is_finite(mov)){ //if1
      
      imod_R(i_n_mod_r)=-1;
      iw=find(imod_R>-1); imod_R=imod_R.elem(iw);
      M_log_RSMP=M_log_RSMP.elem(iw);
          
      limodR=limodR-1;
      t=t-1;
        
     } else{
        
       vtmp1=(M_q+C)/(1-M_q+C);
       vtmp2=(2*(1-G))-1;
       qh=pow_vec(&vtmp1, &vtmp2);
       qh=qh.elem(conv_to<uvec>::from(mov));
    
        
       if(mov.n_elem==1){ //if2
           
        imod_R(i_n_mod_r)=-1;
        iw=find(imod_R>-1); imod_R=imod_R.elem(iw);
        M_log_RSMP=M_log_RSMP.elem(iw);
      
         limodR=limodR-1;  
         edge=mov(0);
         
        } else{
           
           qh.max(imaxe);
           edge=mov(imaxe);
          
        } // end if2
    
        
        G(edge)=1-G(edge);
        add=G(edge);   
       
        
        log_FBF1=log_FBF_Ga_Gb(&G,&G_t,edge,edges,YtY,add,nt,h);
        log_FBF_G=log_FBF1+log_FBF_t;
      
        M_G.col(t)=G;
        
        treeRes=add_to_tree(&G,lM,t,&tree,ltree);
        tree=treeRes(0,0);
        ltree=conv_to<double>::from(treeRes(1,0));
      
        M_log_FBF(t)=log_FBF_G;
        log_pi_G=-log(lcv+1)-lchoose(lcv,sum(G));
        log_num_MP_G=log_FBF_G+log_pi_G;
        log_num_MP(t)=log_num_MP_G;
    
        sum_log_FBF=log_add(sum_log_FBF, log_num_MP_G);
        
        for(imq=0; imq<lcv; imq++){
         if(G(imq)==1){log_sume(imq)=log_add(log_sume(imq), log_num_MP_G);}
        }
        M_q=exp(log_sume-sum_log_FBF);
	    
        limodR=limodR+1;
        imod_R(limodR)=t;
        M_log_RSMP(limodR)=log_num_MP_G; 
         
     } //end if1
    
     t=t+1;
     s=s+1;
    
    } //end while1 

  }// end if univariate 

  t=t-1; 
  s=s-1;
  
  M_P.subvec(0,t)=exp(log_num_MP.subvec(0,t)-sum_log_FBF);
  if(max(M_P.subvec(0,t))>0){
   M_P.subvec(0,t)=M_P.subvec(0,t)/sum(M_P.subvec(0,t));
  }else{
   M_P.subvec(0,t)=zeros<vec>(t+1);	  	  
  }
  
  M_G=M_G.submat(0,0,lcv-1,t);
  
  Res(0,0)=M_q;
  Res(1,0)=M_G;
  Res(2,0)=M_P;
  Res(3,0)=M_log_FBF.subvec(0,t);

  return Res;
   

} // end function



field<mat> FBF_RS(Mat<double>* Corr_c, double nobs_c, Col<double>* G_base_c, double h_c, double C_c, double n_tot_mod_c, double n_hpp_c, bool univariate)
{

 uli neq, rr; 
 double maxne, Mlogbin_sum, lcv, rrmax, q;
 vec V1, V2, vlcv, vG_base, M_q, M_P, iM_P, M_P2;
 mat edges, G_fin, M_G, M_G2;
 mat YtY;
 field<mat> heartRes;
 
 q=(*Corr_c).n_cols; 

 maxne=nobs_c-2*h_c-2;

 neq=1;

 V1=linspace<vec>(1,q-1,q-1);
 V2=zeros<vec>(q-neq);

 edges=join_rows(V1,V2); 
   
 lcv=V1.n_elem;
 vlcv=linspace<vec>(0,lcv-1,lcv); 
  
 vG_base=flipud(*G_base_c);
   
 rrmax=std::min(maxne,lcv); 
 Mlogbin_sum=0;
   
 for(rr=1; rr<=rrmax; rr++){
  Mlogbin_sum=log_add(Mlogbin_sum,lchoose(lcv,rr));   
 }
 
 n_tot_mod_c=std::min(Mlogbin_sum,log(n_tot_mod_c));
 n_tot_mod_c=round(exp(n_tot_mod_c)); 

 YtY=(*Corr_c)*nobs_c;
 heartRes=FBF_heart(nobs_c, &YtY, &vG_base, lcv, &vlcv, &edges, n_tot_mod_c, C_c, maxne, h_c, univariate);
 
 return(heartRes);

} // end function


void printA(string msg)
{

 #ifdef language_cpp
  cout << msg << endl;
 #endif
 
 #ifdef language_py
  msg="print('" + msg + "')";
  PyRun_SimpleString(msg.c_str());
 #endif

 #ifdef language_R	
  Rcout << msg << endl;
 #endif 
	
}


void xit(){
  vector<int> v;
  printA("execution intentionally interrupted");
  printA(to_string(v[0]));
}

template <typename T>
void printV(vector<T> vec,string name){
  printA(name+": ");
  for (auto i: vec){
    printA(to_string(i));
  }
  printA("");
}

template<typename T>
string vec_to_string (T v, uli len){
    
    // string type0=typeid(v(0)).name();
    // bool flg_string=0;
    // if(type0.find("string")!=string::npos){
    //   flg_string=1;
    // }
    
    string res;
    if(len>0){
     res=to_string(v(0));
     for(uli t=1; t<len; ++t){
       res=res+","+to_string(v(t));
     }
    }else{
     res="empty";
    }
    return(res);
}


template <typename T>
string NumberToString ( T Number )
{
   ostringstream ss;
   ss << Number;
   return ss.str();
}

vector<long int> split_string(const string &s, uli order) {
    
	char delim=' ';
	vector<long int> result(order,-1);
    stringstream ss (s);
    string item;

	uli h=0;
    while (getline (ss, item, delim)) {
		result[h]=stoi(item);
		h=h+1;
    }
		
    return result;
}

template<typename T>
uli find_consecutive_finite(T* x, uli col){
   
   uli max_num=0;
   uli num=0;
   for(uli j=0; j<(*x).n_rows; ++j){
    if(isfinite((*x)(j,col))==1){
     num=num+1;
     if(num>max_num){
       max_num=num;
     }
    }else{
     num=0; 
    } 
   }

  return(max_num);

}

template<typename T>
uli find_consecutive_nan(T* x, uli col){
   
   uli max_num=0;
   uli num=0;
   for(uli j=0; j<(*x).n_rows; ++j){
    if(isfinite((*x)(j,col))==0){
     num=num+1;
     if(num>max_num){
       max_num=num;
     }
    }else{
     num=0; 
    } 
   }

  return(max_num);

}


// template<typename T>
// void save_mat(T* obj,string dir, string name, string type){
//   struct timeval time_now{};
//   gettimeofday(&time_now, nullptr); 
//   time_t msecs_time = (time_now.tv_sec * 1000) + (time_now.tv_usec / 1000);
//   (*obj).save(dir + name + "_" + to_string(msecs_time) + "." + type, csv_ascii);
// }


string f_print_perc(double num){ 
 
 string res;
 if(num>=1){
  res=to_string((double)(floor(num*10000)/100)).substr(0,6);    
 }else if(num>=0.1){ 
  res=to_string((double)(floor(num*10000)/100)).substr(0,5); 
 }else{
  res=to_string((double)(floor(num*10000)/100)).substr(0,4);    	   
 } 
 return(res);
}


vector<string> subvector(vector<string> v, Col<uli> idx){
 vector<string> sub_v;
 for(uli j=0; j<idx.n_rows; ++j){
  sub_v.push_back(v[idx(j)]);
 }
 return(sub_v);
}


int rand11(){
 std::random_device rd; 
 std::mt19937 gen(rd()); 
 std::uniform_int_distribution<> distrib(0,  RAND_MAX);
 return(distrib(gen));
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FUNCTIONS FOR VARIABLE SELECTION AND PREDICTIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct str_input
{

  mat MY0;
  mat MY;
  mat corY;
  uli corY_nr;
  uli corY_nc;

  // mat* prt_MY0(){
  //  return(&MY0);
  // }
  
  // mat* prt_MY(){
  //  return(&MY);
  // }

};


str_input data_preparation(mat* Y, Col<uli> vretard)
{

  //varibile target in differenze ritardata
  
  mat MY;
  if(vretard.n_rows>0){
    MY.resize((*Y).n_rows,vretard.n_rows);
      
    for(uli k=0; k<vretard.n_rows; ++k){
     MY.col(k)=shift((*Y).col(0),vretard(k));
     if(vretard(k)>0){
       MY.submat(0,k,vretard(k)-1,k)+=datum::nan;
     }
    }
    //join
    MY=join_rows((*Y),MY);
  }else{
   MY=(*Y);
  }

  //output

  str_input tab_input;
  
  tab_input.MY0=MY;

  MY=MY.rows(find_finite(sum(MY,1)));

  tab_input.MY=MY;

  mat corY=cor(MY);

  tab_input.corY=corY;
  tab_input.corY_nr=MY.n_rows;
  tab_input.corY_nc=MY.n_cols;


  return(tab_input);

}



struct str_output_reg
{
    Col<double> vactual;
    Col<double> vfitted0;
    Col<double> vresid0;
    Col<double> vfitted;
    Col<double> vresid;
    Col<double> vbeta;
    double TSS;
    double RSS;
    double L;
    double L_adj;
};

double f_loss_function(const vec& vals_inp, vec* grad_out, void* opt_data)
{
      
  str_input* data = reinterpret_cast<str_input*>(opt_data);
  mat MY=data->MY;

  double err=0;
  double fit=0;

  uli n=MY.n_rows;

  for(uli t=0; t<n; ++t){
    fit=0;
    for(uli k=0; k<vals_inp.n_rows; ++k){
      fit=fit+vals_inp(k)*MY(t,k+1);
    }
    err=err+(abs(MY(t,0)-fit)/n);
  }

  return err;

}

str_output_reg reg(str_input* tab_input, string loss_function)
{
     
    vec vbeta;

    mat X=(*tab_input).MY.cols(1,(*tab_input).MY.n_cols-1);
    vec Y=(*tab_input).MY.col(0);
    mat XtX=trans(X)*X;
    vec Xty=trans(X)*Y;
    vbeta=solve_linear(&XtX,&Xty); 

    if(loss_function=="MAE"){
     optim::de(vbeta,f_loss_function,tab_input);
    }
    //vec vbeta=inv_sympd(XtX)*Xty;
	    
    vec vactual=Y;
    vec vfitted=X*vbeta;
    vec vresid=vactual-vfitted;
	
    vec vactual0=(*tab_input).MY0.col(0);
    vec vfitted0=(*tab_input).MY0.cols(1,(*tab_input).MY0.n_cols-1)*vbeta;
    vec vresid0=vactual0-vfitted0;
     
    str_output_reg Res;

    Res.vbeta=vbeta;
    Res.vactual=Y;
    Res.vfitted0=vfitted0;
    Res.vresid0=vresid0;
    Res.vfitted=vfitted;
    Res.vresid=vresid;

    return(Res);
}


struct str_out_uni_select
{
 Col<uli> vars_x_idx;
 Col<uli> vars_ar_idx; 
 vec logs_FBF_x;
 vec logs_FBF_ar;
 uli odiff;
 
};

str_out_uni_select model_univariate_selection(mat* Y, double from_lag, double max_lag, bool flg_diff){

  str_out_uni_select str_out;
  Col<uli> vretard;
  uli odiff=0;
  
  if((max_lag!=0) || (((*Y).n_cols-1)>0)){
    
    double cons_rows;
    uli max_lag0;

    cons_rows=(double) find_consecutive_finite(Y,0);

    // if(log(cons_rows/10)>0){
    //  max_lag0=(uli) min(((cons_rows/2)+1),cons_rows/log(cons_rows/10));
    // }else{
    //  max_lag0=1; 
    // }
    if((cons_rows-8)>0){
     max_lag0=cons_rows-8;
    }else{
     max_lag0=1; 
    }

    if(max_lag==-1){
     max_lag=(double) max_lag0;
    }else if(max_lag>0){
     max_lag=(double) min(max_lag,(double)max_lag0);
    }
  
    Col<uli> ids_vars_x;
    vec log_FBF_x;
    Col<uli> ids_vars_ar;
    vec log_FBF_ar;
    
    if(from_lag<=max_lag){
    
      if(max_lag>0){
       vretard=linspace<Col<uli>>((uli)from_lag,(uli)max_lag,(uli)(max_lag-from_lag+1));
      }
  
      double nx=(*Y).n_cols-1;
      double nretard=vretard.n_rows;
      double nvars=nx+nretard;
      
      double h_c=1; 
      if(nvars>50){
       h_c=2; 
      } 
      double n_tot_mod_c=nvars*100;
      double C_c=0.01;
      double threshold; 
      
      str_input tab_input;
      field<mat> res;
    
      uli nr,nc,nvar_max;
      Col<double> G_base_c;
      bool univariate=1;
      Col<double> M_q;
      
      vec v_log_FBF;
      mat M_log_FBF(1,3),M_log_FBF_x,M_log_FBF_ar;
  
      uvec idfinite=find_finite((*Y).col(0));
      nr=idfinite.n_rows;
      
      double qt;
      if(nvars<10){
        qt=0.90;
      }else{
        qt=0.99;
      }
  
  
      vec x;
      vec G_a(1);G_a(0)=1;
      vec G_b(1);G_b(0)=0;
      uli edge=0;
      mat edges(1,2);edges(0,0)=0;edges(0,1)=1;
      uli add=1;
      double log_FBF;
    
      //estimate threshold
      
      uli nv1=(uli)((double)nx);
      uli nv=std::max((uli)200,nv1);
      
      vec y=(*Y).col(0);
      y=y.elem(idfinite);
  
      mat YtY(2,2);YtY(0,0)=1;YtY(1,1)=1;
      v_log_FBF.set_size(nv);
      
      //arma_rng::set_seed(1234567); 
      mt19937 gen(1234567);
      normal_distribution<double> distribution(0.0,1.0);
      x.resize(nr);
      
      for(uli j=0; j<nv; ++j){
        
        //x = randu<vec>(nr,1);
        for(uli k=0; k<nr; ++k){
          x(k)=distribution(gen);
        }
        
        YtY(1,0)=as_scalar(cor(y,x));
        YtY(0,1)=YtY(1,0);
        
        v_log_FBF(j)=log_FBF_Ga_Gb(&G_a, &G_b, edge, &edges, &YtY, add, nr, h_c);
        
      }
      
      uli nth=(uli) (qt*v_log_FBF.size());
      v_log_FBF=sort(v_log_FBF);
      threshold=v_log_FBF(nth);
      
      mat MY((*Y).col(0).n_rows,2);
      MY.col(0)=(*Y).col(0);
      mat MY1;
      

      //select X variables
      
      if(nx>0){
      
        M_log_FBF_x.set_size(nx,3);
        
        for(uli j=0; j<nx; ++j){
          
          MY.col(1)=(*Y).col(j+1);
          MY1=MY.rows(find_finite(sum(MY,1)));
          nr=MY1.n_rows;
          YtY=cor(MY1);
          
          M_log_FBF_x(j,0)=0;
          M_log_FBF_x(j,1)=j+1;
          M_log_FBF_x(j,2)=log_FBF_Ga_Gb(&G_a, &G_b, edge, &edges, &YtY, add, nr, h_c);
          
        }
        
        
        M_log_FBF=join_cols(M_log_FBF,M_log_FBF_x);
      
      }else{
        
        //differentiation
      
        if(flg_diff==1){
        
          double fbf=datum::inf;
          mat Y0;
          
          while(((fbf-threshold)>0) & (odiff<2)){
            MY.col(0)=(*Y).col(0);
            MY.col(1)=linspace<vec>(1,(*Y).n_rows,(*Y).n_rows);
            MY1=MY.rows(find_finite(sum(MY,1)));
            nr=MY1.n_rows;
            YtY=cor(MY1);
            fbf=log_FBF_Ga_Gb(&G_a, &G_b, edge, &edges, &YtY, add, nr, h_c);

            if((fbf-threshold)>0){
              Y0=(*Y);
              (*Y).col(0)=shift((*Y).col(0),1);
              (*Y).submat(0,0,0,0)+=datum::nan;
              for(uli j=0; j<(*Y).n_rows; ++j){
               (*Y)(j,0)=Y0(j,0)-(*Y)(j,0);
              }
              odiff=odiff+1;
            }else{
              break;
            }
          
          }
          
        }
        
      
      }
      
      //select retard
      
      if(nretard>0){
        
        M_log_FBF_ar.set_size(nretard,3);
        
        for(uli j=0; j<nretard; ++j){
          
          MY.col(1)=shift((*Y).col(0),vretard(j));
          if(vretard(j)>0){
            MY.submat(0,1,vretard(j)-1,1)+=datum::nan;
          }
          
          MY1=MY.rows(find_finite(sum(MY,1)));
          nr=MY1.n_rows;
          YtY=cor(MY1);
          
          M_log_FBF_ar(j,0)=1;
          M_log_FBF_ar(j,1)=vretard(j);
          M_log_FBF_ar(j,2)=log_FBF_Ga_Gb(&G_a, &G_b, edge, &edges, &YtY, add, nr, h_c);
        
        }
        
        M_log_FBF=join_cols(M_log_FBF,M_log_FBF_ar);
        
      }
      
      M_log_FBF.shed_row(0);
      
      uvec ids=find(M_log_FBF.col(2)<=threshold);
      M_log_FBF.shed_rows(ids);
      
      
      if(M_log_FBF.n_rows>0){
      
        nvar_max=(uli)2*std::pow(idfinite.n_rows,0.25); //nb: only non missing target is considered
        
        if(M_log_FBF.n_rows>nvar_max){
          ids = sort_index(M_log_FBF.col(2),"descend");
          ids=ids.rows(0,nvar_max-1);
          ids=sort(ids);
          M_log_FBF=M_log_FBF.rows(ids);
        }
      
        uvec ids_x=find(M_log_FBF.col(0)==0);
  
        if(ids_x.size()>0){
          ids_vars_x=conv_to< Col<uli> >::from(M_log_FBF.col(1));
          ids_vars_x=ids_vars_x.elem(ids_x);
          log_FBF_x=M_log_FBF.col(2);
          log_FBF_x=log_FBF_x.elem(ids_x);
        }
      
        uvec ids_ar=find(M_log_FBF.col(0)==1);
      
        if(ids_ar.size()>0){
          ids_vars_ar=conv_to< Col<uli> >::from(M_log_FBF.col(1));
          ids_vars_ar=ids_vars_ar.elem(ids_ar);
          log_FBF_ar=M_log_FBF.col(2);
          log_FBF_ar=log_FBF_ar.elem(ids_ar);
        }
      
      }
    
    }
    
    str_out.vars_x_idx=ids_vars_x; //nb. first regressor has id=1
    str_out.vars_ar_idx=ids_vars_ar;
    str_out.logs_FBF_x=log_FBF_x;
    str_out.logs_FBF_ar=log_FBF_ar;
    str_out.odiff=odiff;
  
  } 
   
  
  return(str_out);

}


field<vec> model_multivariate_selection(mat* Y, Col<uli>* ids_vars_x_uni, Col<uli>* ids_vars_ar_uni, string loss_function, bool flg_const){
        
    Col<uli> ids_vars_x,v_ar;
    vec vbeta_arx;
    vec vresid((*Y).n_rows,1);

    double nvars=(*Y).n_cols-1+(*ids_vars_ar_uni).n_rows; 
    
    if(nvars>0){

      double threshold=0.5;
      double h_c=1; 
      if(nvars>50){
       h_c=2; 
      } 
      double n_tot_mod_c=nvars*100;
      double C_c=0.01;
      double n_hpp_c=1;
       
      mat Y0;
      uvec u_ids_vars_x;
      if((*ids_vars_x_uni).n_rows>0){
       u_ids_vars_x=join_cols(zeros<uvec>(1),conv_to< uvec >::from(*ids_vars_x_uni));
       Y0=(*Y).cols(u_ids_vars_x);
      }else{
       Y0=(*Y).cols(0,0);
      }
    
      str_input tab_input=data_preparation(&Y0,(*ids_vars_ar_uni));
    
      uli nvar_max=(uli)std::pow(tab_input.MY.n_rows,0.25);
      
      uli nr,nc;
      Col<double> G_base_c;
      bool univariate=0;
      vec M_q;
      field<mat> res;
      uvec ids;
      vec rvals;
  
      uli len_ar=(*ids_vars_ar_uni).n_rows;
      vec ids_vars;
  
      mat M_G;
      vec M_P;
    
      if((len_ar>0) || (Y0.n_cols>1)){
    	
	      nr=tab_input.corY_nr;
        nc=tab_input.corY_nc;
  
        h_c=1; 
        if((nc-1)>50){
         h_c=2; 
        } 
    
        G_base_c=zeros<vec>(nc-1);
        
        res=FBF_RS(&tab_input.corY,nr,&G_base_c,h_c,C_c,n_tot_mod_c,n_hpp_c,univariate);
      
                
        M_q=res(0,0);
  
	      ids=find(M_q>=threshold);
	      if(ids.n_rows>nvar_max){
          ids = sort_index(M_q,"descend");
          ids=ids.rows(0,nvar_max-1);
          ids=sort(ids);
        }
    
        rvals=linspace<vec>(0,M_q.n_rows-1, M_q.n_rows).elem(ids);
      
        ids_vars=rvals-(Y0.n_cols-1);
        vec ids_vars_ar=ids_vars(find(ids_vars>=0));
        uvec fd_ids_vars=find(ids_vars<0);
        if(fd_ids_vars.n_cols>0){
         ids_vars_x=conv_to< Col<uli> >::from(rvals(fd_ids_vars));
        }
    
	      if(ids_vars.n_rows>0){
          if(len_ar>0){
            uvec u_ids_vars_ar=conv_to< uvec >::from(ids_vars_ar);
            v_ar=(*ids_vars_ar_uni).rows(u_ids_vars_ar);
          }
          
          uvec urvals=join_cols(zeros<uvec>(1),conv_to< uvec >::from(rvals+1));
          tab_input.MY0=tab_input.MY0.cols(urvals);
          tab_input.MY=tab_input.MY.cols(urvals);
        
          if(flg_const==1){
           tab_input.MY0.resize(tab_input.MY0.n_rows,tab_input.MY0.n_cols+1);
           tab_input.MY0.col(tab_input.MY0.n_cols-1)=zeros<vec>(tab_input.MY0.n_rows)+1;

           tab_input.MY.resize(tab_input.MY.n_rows,tab_input.MY.n_cols+1);
           tab_input.MY.col(tab_input.MY.n_cols-1)=zeros<vec>(tab_input.MY.n_rows)+1;
          }
          
          str_output_reg tab_out_reg=reg(&tab_input, loss_function);
          vbeta_arx=tab_out_reg.vbeta;
          vresid=tab_out_reg.vresid0;
        }else{
	    	  vresid=Y0.col(0);
          double m0=as_scalar(mean(vresid.rows(find_finite(vresid))));
          for(uli t=0; t<vresid.n_rows; ++t){
            vresid(t)=vresid(t)-m0;
	        }
        }
      }else{
    
        vresid=Y0.col(0);
        double m0=as_scalar(mean(vresid.rows(find_finite(vresid))));
        for(uli t=0; t<vresid.n_rows; ++t){
         vresid(t)=vresid(t)-m0;
        }
      
      }

      if(ids_vars_x.size()>0){ 
       u_ids_vars_x=conv_to< uvec >::from(ids_vars_x);
       ids_vars_x=(*ids_vars_x_uni).elem(u_ids_vars_x);
      }
      
    }else{

     vresid=(*Y).col(0);

    }

    
    uli nmodels=1;
    
    field<vec> str_out;
    str_out.set_size(nmodels,4);
    str_out(0,0)=conv_to< vec >::from(ids_vars_x);
    str_out(0,1)=conv_to< vec >::from(v_ar);
    str_out(0,2)=vbeta_arx;
    str_out(0,3)=vresid;
    
    return(str_out);

}



map<string,double> performances(vec vactual, vec vfitted, uli nvars){

  vec vresid=vactual-vfitted;
  uvec non_missing=find_finite(vresid);
  vresid=vresid.elem(non_missing);   
  //double RSS_=as_scalar(sum(pow(vresid,2)));
  
  vactual=vactual.elem(non_missing);
  double m0=mean(vactual);
  // Col<double> vTSS=vactual;
  // for(uli i=0; i<vTSS.n_rows; ++i){
  //  vTSS(i)=vTSS(i)-m0;
  // }
  // vTSS=pow(vTSS,2);
  // double TSS=sum(vTSS);

  // double R2 = 1 - (RSS/TSS);
  // double R2_adj = 1 - (((double)vactual.n_rows-1)/((double)vactual.n_rows-(double)nvars-1))*(RSS/TSS); 
 
  //abs dist


  double L1=as_scalar(sum(abs(vresid)));

  double L0=0;
  for(uli i=0; i<vactual.n_rows; ++i){
   L0=L0+abs(vactual(i)-m0);
  }  

  double L=1-(L1/L0);

  double L_adj=1 - (((double)vactual.n_rows-1)/((double)vactual.n_rows-(double)nvars-1))*(L1/L0);

  map<string,double> res;
  res["L"]=L;
  res["L_adj"]=L_adj;

  return(res);
    
}


struct str_pred_out
{

 mat predictions;
 mat fitted;
 double L;
 double L_adj;
 
};


str_pred_out sarimax_pred(mat* Y, bool flg_sim, mat Mfitted, vec probs, uli nsim, string loss_function, bool pred_only, bool flg_const, field<vec> models)
{
  
  Col<uli> ids_vars_x=conv_to< Col<uli> >::from(models(0,0));
  Col<uli> ids_vars_ar=conv_to< Col<uli> >::from(models(0,1));
  vec vbeta_arx=models(0,2);
  ///
  double const0=0;
  if(flg_const==1){
   if(vbeta_arx.n_rows>0){
    const0=vbeta_arx(vbeta_arx.n_rows-1);
    vbeta_arx.shed_row(vbeta_arx.n_rows-1);
   }
  }
  ///
  
  Col<uli> ids_vars_ma=conv_to< Col<uli> >::from(models(0,3));
  vec vbeta_ma=models(0,4);
  
  str_pred_out str_out;
  vec vresid; 

  uli p=ids_vars_ar.n_rows;
  uli q=ids_vars_ma.n_rows;
  uli k=vbeta_arx.n_rows-p; //number of regressors

  uli maxpq=0;
  if((p==0) && (q!=0)){
    maxpq=ids_vars_ma.max();
  }else if((p!=0) && (q==0)){
    maxpq=ids_vars_ar.max();
  }else if((p!=0) && (q!=0)){
    maxpq=max(ids_vars_ar.max(),ids_vars_ma.max());
  }
  
  uli npred;
  if(flg_sim==0){
    uli npred1=find_consecutive_nan(Y,0);
    uli npred2=find_consecutive_finite(Y,0)-maxpq-k;

    npred=min(npred1,npred2);
    if(npred==0){
      npred=1;
    }
  }else{
    npred=Mfitted.n_cols;
  }
  
  uli ri;
  double rd;
  // //random_device rdv; 
  mt19937 gen(1234567);
  // default_random_engine gen;
  uniform_real_distribution<double> distrib(0.0,1.0);

  // mt19937 engine(1234567);
  // uniform_real_distribution<double> distrib(0.0,1.0);
  // auto gen = bind(ref(distrib), ref(engine));

  field<vec> Fresid(npred,1);
  vec resid_size(npred);

  if(flg_sim==1){ 
     for(uli kpred=0; kpred<npred; ++kpred){
       vresid=(*Y).col(0)-Mfitted.col(kpred);
       vresid=vresid(find_finite(vresid));
       Fresid(kpred,0)=vresid;
       resid_size(kpred)=(double) vresid.size();
     }
  }else{
    nsim=1;
  }

  vec vbeta_x;
  vec vbeta_ar;

  if(k>0){
   vbeta_x=vbeta_arx.rows(0,k-1);
  }
  if(p>0){
   vbeta_ar=vbeta_arx.rows(k,vbeta_arx.n_rows-1);
  }
  double tar;
  double tma;
  
  mat My_pred((*Y).n_rows, npred);
  vec veps((*Y).n_rows);
  field<mat> Fout(npred,1);
  mat Mout;
  
  if(flg_sim==1){
    for(uli kpred=0; kpred<npred; ++kpred){
     Fout(kpred,0).resize(nsim,(*Y).n_rows);
     Fout(kpred,0).fill(datum::nan);
    }
    }else{  
     for(uli kpred=0; kpred<npred; ++kpred){
      Fout(kpred,0).resize((*Y).n_rows,2);
    }
  }
  
  double pred_x=0;
  double pred_ar=0;
  double pred_ma=0;
  double pred_err;
  
  uvec ut(1);
  uvec uids_vars_x=conv_to<uvec>::from(ids_vars_x);

  double eps_tma;

  bool is_na_pred_ar=0;
  
  uli t0, tk, last_kpred;

  for(uli s=0; s<nsim; ++s){

    My_pred.fill(datum::nan);
    veps.fill(datum::nan);

    if(maxpq>0){
     t0=(maxpq+1);
    }else{
     t0=0; 
    }
      

    for(uli t=t0; t<(*Y).n_rows; ++t){
      
      last_kpred=std::min(npred,(uli)(*Y).n_rows-t);
      
      for(uli kpred=0; kpred<last_kpred; ++kpred){ 

        tk=t+kpred;

        ut(0)=tk;
   
        //X
        
        if(k>0){
         pred_x=as_scalar((*Y).submat(ut,uids_vars_x)*vbeta_x);
        }else{
         pred_x=0;
        }
        
        //AR

        pred_ar=const0;
    
        if(p>0){
        
          for(uli p0=0; p0<p; ++p0){
           
           tar=(double)tk-(double)ids_vars_ar(p0);
           
           if(tar>0){
            if((tar<t) && (isfinite((*Y)(tar,0)))){ 
             pred_ar=pred_ar+(*Y)(tar,0)*vbeta_ar(p0);
            }else{ 
             is_na_pred_ar=1;
             for(uli z=0; z<=kpred; ++z){ 
              if(isfinite(My_pred(tar,z))){ 
               pred_ar=pred_ar+My_pred(tar,z)*vbeta_ar(p0);
               is_na_pred_ar=0;
               break;
              }
             }
             if(is_na_pred_ar==1){
               pred_ar=datum::nan;
             }
            }
           }else{
            pred_ar=datum::nan;
            break; 
           }
      
          }
        
        }//end if p
    
        pred_ma=0;

        if(isfinite(pred_ar)){
        
          //MA
      
          if(q>0){
              
            for(uli q0=0; q0<q; ++q0){
        
             tma=(double)tk-(double)ids_vars_ma(q0);
        
             if(tma>0){
              eps_tma=veps(tma);
              if((tma<t) && (isfinite(eps_tma))){
                pred_ma=pred_ma+eps_tma*vbeta_ma(q0);
              }else{
                pred_ma=pred_ma+0;
              }
             }else{
              pred_ma=pred_ma+0;
             }
        
            }
            
          }//end if q

    
        }
            
        if(isfinite(pred_ar)){
         My_pred(tk,kpred)=pred_x+pred_ar+pred_ma; 
        }else{
         My_pred(tk,kpred)=datum::nan; 
        }
            
        if(flg_sim==1){
                    
          rd=distrib(gen);
          ri=(uli)floor(resid_size(kpred)*rd);

          pred_err=Fresid(kpred,0)(ri);

          if(kpred==0){
            if(isfinite((*Y)(t,0))){
              veps(t)=(*Y)(t,0)-pred_x-pred_ar;
            }else{
              veps(t)=0+pred_err;
            }
          }
          
          if(isfinite(pred_ar)){
           My_pred(tk,kpred)=pred_x+pred_ar+pred_ma+pred_err; 
          }else{
           My_pred(tk,kpred)=datum::nan; 
          }

          Fout(kpred,0)(s,tk)=My_pred(tk,kpred);

        }else{ 

          if(kpred==0){
           if(isfinite((*Y)(t,0))){
            veps(t)=(*Y)(t,0)-pred_x-pred_ar;
           }else{
            veps(t)=0;
           }
          }
        
        }
        
        // if(kpred==0){
        //  cout << "--------" << endl;
        //  cout << "tk:" << tk << endl;
        //  cout << "pred_ma:" << pred_ma << endl;
        //  cout << "pred_ar:" << pred_ar << endl;
        //  cout << "eps-1:" << veps(tk-1) << endl;
        //  cout << "eps:" << veps(tk) << endl;
        //  cout << "target-1:" << (*Y)(tk-1,0) << endl;
        //  cout << "target:" << (*Y)(tk,0) << endl;
        //  cout << "pred:" << My_pred(tk,0) << endl;
        // }
      
      }//end kpred 
      
    }//end t
  
  }//end s
  
  if(flg_sim==1){
   
   for(uli kpred=0; kpred<npred; ++kpred){ 
     Fout(kpred,0)=quantile(Fout(kpred,0),probs);
     Fout(kpred,0)=join_rows(Mfitted.col(kpred),Fout(kpred,0).t());
     Fout(kpred,0)=join_rows((*Y).col(0),Fout(kpred,0));
   }

   if(npred>1){
      
     //prediction interval correction

     uli kpred=0;
     uli h_low=0, h_up=0;
   
     mat M_low((*Y).n_rows,3); M_low.fill(datum::nan);
     mat M_up((*Y).n_rows,3); M_up.fill(datum::nan);

     mat M0_low((*Y).n_rows,3); M0_low.fill(datum::nan);
     M0_low.col(1)=Fout(0,0).col(2);
     M0_low.col(2).fill(1);
     mat M0_up((*Y).n_rows,3); M0_up.fill(datum::nan);
     M0_up.col(1)=Fout(0,0).col(4);
     M0_up.col(2).fill(1);

     for(uli j=0; j<(*Y).n_rows; ++j){
       if(isfinite((*Y)(j,0))==0){
         if((kpred>0) && (kpred<npred)){
           M_low(h_low,0)=Fout(kpred,0)(j,2);
           M_low(h_low,1)=Fout(0,0)(j,2);
           M_low(h_low,2)=1;
           M0_low(j,0)=Fout(kpred,0)(j,2);
           h_low=h_low+1;
  
           M_up(h_up,0)=Fout(kpred,0)(j,4);
           M_up(h_up,1)=Fout(0,0)(j,4);
           M_up(h_up,2)=1;
           M0_up(j,0)=Fout(kpred,0)(j,4);
           h_up=h_up+1;
         }
         kpred=kpred+1;  
       }else{
        kpred=0; 
       }  
     }

     str_input tab_input_ci;
     str_output_reg reg_ci; 
     vec vlower_bound, vupper_bound;
     vec vquantiles_ci;

     tab_input_ci.MY0=M0_low;
     tab_input_ci.MY=M_low.rows(0,h_low-1);
     reg_ci=reg(&tab_input_ci, loss_function);
     vquantiles_ci=quantile(reg_ci.vresid,probs);
     vlower_bound=reg_ci.vfitted0;
     for(uli j=0; j<vlower_bound.n_rows; ++j){
      vlower_bound(j)=vlower_bound(j)+vquantiles_ci(0);
     }

     tab_input_ci.MY0=M0_up;
     tab_input_ci.MY=M_up.rows(0,h_low-1);
     reg_ci=reg(&tab_input_ci, loss_function);
     vquantiles_ci=quantile(reg_ci.vresid,probs);
     vupper_bound=reg_ci.vfitted0;
     for(uli j=0; j<vupper_bound.n_rows; ++j){
      vupper_bound(j)=vupper_bound(j)+vquantiles_ci(2);
     }  

     Mout=Fout(0,0);

     kpred=0;
     for(uli j=0; j<(*Y).n_rows; ++j){
       if(isfinite((*Y)(j,0))==0){
         if((kpred>0) && (kpred<npred)){
           if((isfinite(vlower_bound(j))) && isfinite(vupper_bound(j))){
            Mout(j,2)=vlower_bound(j);
            Mout(j,4)=vupper_bound(j);
           }
         }
         kpred=kpred+1;  
       }else{
        kpred=0; 
       }  
     }

   }else{

    Mout=Fout(0,0);

   }//end if(npred>1)

  }else{
  
   uvec nonmiss;

   for(uli kpred=0; kpred<npred; ++kpred){
     nonmiss=find_finite((*Y).col(0)-My_pred.col(kpred));
     npred=kpred+1;
     if(nonmiss.n_rows<20){ //20 is the number of observations over that the forecasting horizon can be considered significative
       break;
     }
   }
  
  }
     
  map<string,double> mp_idx_perf;

  double L, L_adj;
  if(pred_only==0){
   mp_idx_perf=performances((*Y).col(0), My_pred.col(0), (uli)(ids_vars_x.n_rows+ids_vars_ar.n_rows+ids_vars_ma.n_rows));
   L=mp_idx_perf["L"];
   L_adj=mp_idx_perf["L_adj"];
  }
  
  str_out.predictions=Mout;
  if(flg_sim==0){
   str_out.fitted=My_pred.cols(0,npred-1);
  }else{
   str_out.fitted=My_pred; 
  }
  str_out.L=L;
  str_out.L_adj=L_adj;

  return(str_out);

}


struct str_model_out
{
  Col<uli> ids_vars_x;
  Col<uli> ids_vars_ar;
  Col<double> vbeta_arx;
  Col<uli> ids_vars_ma;
  Col<double> vbeta_ma;
};  


struct str_model_selection
{
 
 field<vec> models;
 mat predictions;

};


str_model_selection model_selection_prediction(mat* Y, double from_lag, double max_lag, vec probs, uli nsim, string loss_function, bool pred_only, bool flg_const, bool flg_diff, field<vec> models)
{

  str_pred_out out_pred_i;

  str_out_uni_select out_uni_select_arx;
  field<vec> out_multi_select_arx;

  str_out_uni_select out_uni_select_ma;
  field<vec> out_multi_select_ma;
  
  vec vresid;
  mat Mfitted, Mfitted_empty;
  Col<uli> vretard_empty;

  bool flg_x_only=0;
  if(max_lag==0){
   flg_x_only=1; 
  }
  
  str_model_selection res_out;
  
  res_out.models.set_size(1,6);
  
  //SARIX
  
  uli ndiff=0;
  vec old_target=(*Y).col(0);
  
  if(pred_only==0){
    
    out_uni_select_arx=model_univariate_selection(Y, from_lag, max_lag, flg_diff);
    
    res_out.models(0,5)=out_uni_select_arx.odiff;
    ndiff=(uli)out_uni_select_arx.odiff;
    
  }else{
    
    ndiff=(uli)as_scalar(models(0,5));
    mat Y0;
    for(uli k=0; k<ndiff; ++k){
      Y0=(*Y);
      (*Y).col(0)=shift((*Y).col(0),1);
      (*Y).submat(0,0,0,0)+=datum::nan;
      for(uli j=0; j<(*Y).n_rows; ++j){
       (*Y)(j,0)=Y0(j,0)-(*Y)(j,0);
      }
    }
    
  }
  
  if(((flg_x_only==0) && (out_uni_select_arx.vars_ar_idx.n_rows>0) && (pred_only==0)) || ((flg_x_only==0) && (pred_only==1))){ //if it is a sarimax
    
    
    if(pred_only==0){
      
      
      out_multi_select_arx=model_multivariate_selection(Y, &out_uni_select_arx.vars_x_idx, &out_uni_select_arx.vars_ar_idx, loss_function, flg_const);
      
      
      res_out.models(0,0)=out_multi_select_arx(0,0);
      res_out.models(0,1)=out_multi_select_arx(0,1);
      res_out.models(0,2)=out_multi_select_arx(0,2);
      
      //MA

      vresid=out_multi_select_arx(0,3);

      out_uni_select_ma=model_univariate_selection(&vresid, from_lag, max_lag, 0);
      
      out_multi_select_ma=model_multivariate_selection(&vresid, &out_uni_select_ma.vars_x_idx, &out_uni_select_ma.vars_ar_idx, loss_function, flg_const);
      
      res_out.models(0,3)=out_multi_select_ma(0,1);
      res_out.models(0,4)=out_multi_select_ma(0,2);
      
      models=res_out.models;
    
    }
    
    out_pred_i=sarimax_pred(Y, 0, Mfitted_empty, probs, nsim, loss_function, pred_only, flg_const, models);
    Mfitted=out_pred_i.fitted;

    out_pred_i=sarimax_pred(Y, 1, Mfitted, probs, nsim, loss_function, pred_only, flg_const, models);
    
    res_out.predictions=out_pred_i.predictions;
    
  }else{
    
    flg_x_only=1;
  
  }
  

  if(flg_x_only==1){
     
    if(pred_only==0){
    
      out_multi_select_arx=model_multivariate_selection(Y, &out_uni_select_arx.vars_x_idx, &vretard_empty, loss_function, flg_const);

      res_out.models(0,0)=out_multi_select_arx(0,0);
      res_out.models(0,1)=out_multi_select_arx(0,1);
      res_out.models(0,2)=out_multi_select_arx(0,2);
      
      models=res_out.models;
      
    }
    
    out_pred_i=sarimax_pred(Y, 0, Mfitted_empty, probs, nsim, loss_function, pred_only, flg_const, models);
    Mfitted=out_pred_i.fitted;
    
    out_pred_i=sarimax_pred(Y, 1, Mfitted, probs, nsim, loss_function, pred_only, flg_const, models);
    
    res_out.predictions=out_pred_i.predictions;

  }
  
  if(ndiff>0){
   
      double y_t1,y_t2;
      out_pred_i.predictions.col(0)=old_target;
      
      for(uli j=0; j<ndiff; ++j){
        for(uli k=1; k<5; ++k){
          out_pred_i.predictions(j,k)=datum::nan;
        }
      }
      
      for(uli j=ndiff; j<out_pred_i.predictions.n_rows; ++j){
        for(uli k=1; k<5; ++k){
          if(isfinite(out_pred_i.predictions(j,k))==1){
            y_t1=out_pred_i.predictions(j-1,0);
            if(isfinite(y_t1)==0){
             y_t1=out_pred_i.predictions(j-1,k);
            }
            if(ndiff==1){
             out_pred_i.predictions(j,k)=out_pred_i.predictions(j,k)+y_t1;
            }
            if(ndiff==2){
             y_t2=out_pred_i.predictions(j-2,0);
             if(isfinite(y_t2)==0){
              y_t2=out_pred_i.predictions(j-2,k);
             }
             out_pred_i.predictions(j,k)=out_pred_i.predictions(j,k)+2*y_t1-y_t2;
            }
          }
        }
      }
      
      res_out.predictions=out_pred_i.predictions;
    
  }
  
  (*Y).col(0)=old_target;
  
  return(res_out);

}


struct str_output
{

 mat predictions;
 vec performances;
 
 mat fw_predictions;
 field<vec> fw_models;
 vec fw_performances;
 
 mat bw_predictions;
 field<vec> bw_models;
 vec bw_performances;

};



str_output regpred_cpp(mat* Y, double from_lag, double max_lag, double alpha, uli nsim, bool flg_print, string direction, string loss_function, bool pred_only, bool flg_const, bool flg_diff, vector < field<vec> > vmodels)
{
  
  str_model_out tab_model;
  str_pred_out tab_pred;
  str_output str_out;
   
  double pinf=(alpha/2);
  double psup=1-(alpha/2);

  vec probs={pinf, 0.5, psup};

  vec vresid_empty;
  vec vresid;
  vec vfitted;

  Col<uli> vtmp_uli;
  map<string,double> mp_idx_perf;

  uli nrows=(*Y).n_rows;
  
  mat Yr;
  if((direction=="<->") || (direction=="<-")){ //do not move below
   Yr=reverse((*Y),0);
  }

  uli bw_k=0, fw_k=0;

  mat predictions, predictions_rev;

  str_model_selection out_sel_pred;
  
  field<vec> models0;
  field<vec> models;
  
  double L,L_adj;
  
  if((direction=="<->") || (direction=="->")){ 
    
    //model selection
    if(flg_print==1){
     printA("Forward prediction: model selection and prediction...");
    }
    
    if(pred_only==1){
     models0=vmodels[0];
    }
    
    out_sel_pred=model_selection_prediction(Y, from_lag, max_lag, probs, nsim, loss_function,pred_only,flg_const,flg_diff,models0);
  
    predictions=out_sel_pred.predictions;

    str_out.fw_predictions=predictions;
    
    models=out_sel_pred.models;
    str_out.fw_models=models;
    
    L=datum::nan;
    L_adj=datum::nan;
    if(pred_only==0){
      fw_k=(uli) (models(0,0).size() + models(0,1).size() + models(0,3).size());
      mp_idx_perf=performances(predictions.col(0), predictions.col(3), fw_k);
      L=mp_idx_perf["L"];
      L_adj=mp_idx_perf["L_adj"];
    }
    
    str_out.fw_performances.set_size(2);
    str_out.fw_performances(0)=L;
    str_out.fw_performances(1)=L_adj;
    
    
  }
 
  if((direction=="<->") || (direction=="<-")){ 
                    
    //model selection
    if(flg_print==1){
     printA("Backward prediction: model selection and prediction...");
    }
    
    if(pred_only==1){
     models0=vmodels[1];
    }
    
    out_sel_pred=model_selection_prediction(&Yr, from_lag, max_lag, probs, nsim, loss_function,pred_only,flg_const,flg_diff,models0);

    predictions_rev=reverse(out_sel_pred.predictions);
    
    str_out.bw_predictions=predictions_rev;
    
    models=out_sel_pred.models;
    str_out.bw_models=models;

    L=datum::nan;
    L_adj=datum::nan;
    if(pred_only==0){
      bw_k=(uli) (models(0,0).size() + models(0,1).size() + models(0,3).size());
      mp_idx_perf=performances(predictions_rev.col(0), predictions_rev.col(3), bw_k);
      L=mp_idx_perf["L"];
      L_adj=mp_idx_perf["L_adj"];
    }
    
    str_out.bw_performances.set_size(2);
    str_out.bw_performances(0)=L;
    str_out.bw_performances(1)=L_adj;

  }

  //collapse
  
  if(direction=="<-"){
    predictions=predictions_rev;
  }

  
  if(direction=="<->"){

   for(uli t=0; t<nrows; ++t){
    for(uli k=1; k<4; ++k){
     if(isfinite(predictions(t,1)) && isfinite(predictions_rev(t,1))){
      predictions(t,1)=(predictions(t,1)+predictions_rev(t,1))/2;
     }else if(isfinite(predictions_rev(t,1))){
      predictions(t,1)=predictions_rev(t,1);
     } 
    }
   }

  }

  str_out.predictions=predictions;
         
  L=datum::nan;
  L_adj=datum::nan;
  if(pred_only==0){
    mp_idx_perf=performances(predictions.col(0), predictions.col(3), max(bw_k,fw_k));
    L=mp_idx_perf["L"];
    L_adj=mp_idx_perf["L_adj"];
  }
  
  str_out.performances.set_size(2);
  str_out.performances(0)=L;
  str_out.performances(1)=L_adj;
  
  if(flg_print==1){ 
   printA("Process ended successfully!");
  }
  
  return(str_out);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FUNCTION FOR PASSING RESULTS TO PYTHON AND R 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

mat  std_vec2_to_arma_mat(svec2* A) {

    uli nrows=(*A).size();
    uli ncols=(*A)[0].size();
    
    mat V(nrows,ncols);
    
    for (uli j = 0; j < nrows; ++j) {
     for (uli i = 0; i < ncols; ++i) {
       V(j,i)=(*A)[j][i];
     }
    }
   
    return V;

}

svec2  arma_mat_to_std_vec2(mat* A) {
    
    svec2  V((*A).n_rows);
    for (size_t i = 0; i < (*A).n_rows; ++i) {
        V[i] = conv_to< svec1 >::from((*A).row(i));
    };
    
    return V;
}


field<vec> std_vec3_to_arma_fie_vec(svec3* A){
  
    uli nrows=(*A).size(); //number of models
    uli ncols=(*A)[0].size();; //number of parameter vectors
    
    field<vec> V(nrows,ncols);
    
    for (uli j = 0; j < nrows; ++j) {
     for (uli i = 0; i < ncols; ++i) {
       V(j,i)=(*A)[j][i];
     }
    }
   
    return V;
  
}


svec3 arma_fie_vec_to_std_vec3(field<vec>* A){
  
    uli nrows=(*A).n_rows; //number of models
    uli ncols=(*A).n_cols; //number of parameter vectors
    
    svec3 V(nrows, vector< vector<double> >(ncols , vector<double>()));
    // svec3 V;
    // svec3.resize(nrows);
    
    for (uli j = 0; j < nrows; ++j) {
     //V[j].resize(ncols);
     for (uli i = 0; i < ncols; ++i) {
       //V[j][i].push_back(conv_to< svec1 >::from((*A)(j,i)));
       V[j][i]=conv_to< svec1 >::from((*A)(j,i));
     }
    }
   
    return V;
  
}

#ifdef language_py


pair < svec3, pair < svec2, svec4 > > 
regpred_py(svec2& Y, double from_lag, double max_lag, double alpha, uli nsim, bool flg_print, string direction, string loss_function, bool pred_only, bool flg_const, bool flg_diff, svec4& vmodels){
  
  mat Y0=std_vec2_to_arma_mat(&Y);
  
  vector < field<vec> > vmodels0(2);
  if(pred_only==1){
   if((direction=="->") || (direction=="<->")){
    vmodels0[0]=std_vec3_to_arma_fie_vec(&vmodels[0]);
   }
   if((direction=="<-") || (direction=="<->")){
    vmodels0[1]=std_vec3_to_arma_fie_vec(&vmodels[1]);
   }
  }
    
  str_output str_out=regpred_cpp(&Y0, from_lag, max_lag, alpha, nsim, flg_print, direction, loss_function, pred_only, flg_const, flg_diff, vmodels0);

  //store predictions
    
  svec3 vpredictions(3);
  
  svec2 predictions=arma_mat_to_std_vec2(&str_out.predictions);
  svec2 fw_predictions=arma_mat_to_std_vec2(&str_out.fw_predictions);
  svec2 bw_predictions=arma_mat_to_std_vec2(&str_out.bw_predictions);
  
  vpredictions[0]=predictions;
  vpredictions[1]=fw_predictions;
  vpredictions[2]=bw_predictions;
  
  //store performances
  
  svec2 vperformances(3);
  
  svec1 performances=conv_to< svec1 >::from(str_out.performances);
  svec1 fw_performances=conv_to< svec1 >::from(str_out.fw_performances);
  svec1 bw_performances=conv_to< svec1 >::from(str_out.bw_performances);
  
  vperformances[0]=performances;
  vperformances[1]=fw_performances;
  vperformances[2]=bw_performances;
  
  //store models
  
  svec4 vmodels1(2);
  if(pred_only==0){
     vmodels1[0]=arma_fie_vec_to_std_vec3(&str_out.fw_models);
     vmodels1[1]=arma_fie_vec_to_std_vec3(&str_out.bw_models);
  }
  
  //final store
  
  pair < svec3, pair < svec2, svec4 > > res;
  res.first=vpredictions;
  res.second.first=vperformances;
  res.second.second=vmodels1;
  
  return(res);

}

#endif

#ifdef language_R

template <typename T>
NumericVector arma_vec_to_R_vec(const T* x) {
    return NumericVector((*x).begin(), (*x).end());
}

vec R_vec_to_arma_vec(NumericVector* V) {
    
    uli nr=(uli) (*V).size();
    
    vec O(nr);
    
    for (size_t j = 0; j < nr; ++j) {
      O(j)=(*V)(j);
    }
    
    return(O);
    
}

NumericMatrix  arma_mat_to_R_mat(mat* A) {
    
  NumericMatrix  V((*A).n_rows,(*A).n_cols);
    
	for (size_t j = 0; j < (*A).n_rows; ++j) {
     for (size_t i = 0; i < (*A).n_cols; ++i) {
        V(j,i) =(*A)(j,i);
     }
	}
    
  return V;
}


field<vec> R_List2_vec_to_arma_fie_vec(List L0){
  
  uli nrows=L0.size();
  List L1=L0[0];
  uli ncols=L1.size();
  
  NumericVector v;
  field<vec> O(nrows,ncols);
  
  for (uli i0 = 0; i0 < nrows; ++i0){
    L1=L0[i0];
    for (uli i1 = 0; i1 < ncols; ++i1){
      v=L1[i1];
      O(i0,i1)=R_vec_to_arma_vec(&v);
    }
  }
  
  return(O);

}


List arma_fie_vec_to_R_List2_vec(field<vec>* F){
  
  uli nrows=(*F).n_rows;
  uli ncols=(*F).n_cols;
  
  vec v;
  List R0(nrows); 
  
  for (uli i0 = 0; i0 < nrows; ++i0){
    List R1(ncols);
    for (uli i1 = 0; i1 < ncols; ++i1){
      v=(*F)(i0,i1);
      R1[i1]=arma_vec_to_R_vec(&v);
    }
    R0[i0]=R1;
  }
  
  return(R0);
  
}

RcppExport SEXP regpred_R(SEXP Y_p, SEXP from_lag_p, SEXP max_lag_p, SEXP alpha_p, SEXP nsim_p, SEXP flg_print_p, SEXP direction_p, SEXP loss_function_p, SEXP pred_only_p, SEXP flg_const_p, SEXP flg_diff_p, SEXP vmodels_p)
{

  NumericMatrix Y_0(Y_p); 
  mat Y(Y_0.begin(), Y_0.nrow(), Y_0.ncol(), false);
  
  NumericVector from_lag_0(from_lag_p); 
  double from_lag = Rcpp::as<double>(from_lag_0);
  
  NumericVector max_lag_0(max_lag_p); 
  double max_lag = Rcpp::as<double>(max_lag_0);
  
  NumericVector alpha_0(alpha_p); 
  double alpha = Rcpp::as<double>(alpha_0);
  
  NumericVector nsim_0(nsim_p); 
  uli nsim = Rcpp::as<uli>(nsim_0);
  
  NumericVector flg_print_0(flg_print_p); 
  bool flg_print = Rcpp::as<bool>(flg_print_0);

  CharacterVector direction_0(direction_p); 
  string direction = Rcpp::as<string>(direction_0);

  CharacterVector loss_function_0(loss_function_p); 
  string loss_function = Rcpp::as<string>(loss_function_0);
  
  NumericVector pred_only_0(pred_only_p); 
  bool pred_only = Rcpp::as<bool>(pred_only_0);
  
  NumericVector flg_const_0(flg_const_p);
  bool flg_const = Rcpp::as<bool>(flg_const_0);
  
  NumericVector flg_diff_0(flg_diff_p);
  bool flg_diff = Rcpp::as<bool>(flg_diff_0);
  
  List vmodels(vmodels_p);
  
  vector < field<vec> > vmodels0(2);
  List fw_models;
  List bw_models;
  if(pred_only==1){
   vmodels0[0]=R_List2_vec_to_arma_fie_vec(vmodels[0]);
   vmodels0[1]=R_List2_vec_to_arma_fie_vec(vmodels[1]);
   fw_models=vmodels[0];
   bw_models=vmodels[1];
  }

  str_output str_out=regpred_cpp(&Y, from_lag, max_lag, alpha, nsim, flg_print, direction, loss_function, pred_only, flg_const, flg_diff, vmodels0);
  
  //store predictions
  
  NumericMatrix predictions=arma_mat_to_R_mat(&str_out.predictions);
  NumericMatrix fw_predictions=arma_mat_to_R_mat(&str_out.fw_predictions);
  NumericMatrix bw_predictions=arma_mat_to_R_mat(&str_out.bw_predictions);
  
  //store performances
  
  NumericVector performances=arma_vec_to_R_vec(&str_out.performances);
  NumericVector fw_performances=arma_vec_to_R_vec(&str_out.fw_performances);
  NumericVector bw_performances=arma_vec_to_R_vec(&str_out.bw_performances);
  
  //store models
  
  if(pred_only==0){
     fw_models=arma_fie_vec_to_R_List2_vec(&str_out.fw_models);
     bw_models=arma_fie_vec_to_R_List2_vec(&str_out.bw_models);
  }
  
  /*maximum 20 elements admitted for each level*/
  List res=List::create(
    Named("prediction")=List::create(
      Named("final") = predictions,
      Named("forward") = fw_predictions,
      Named("backward") = bw_predictions
    ),
    Named("performance")=List::create(
      Named("final") = performances,
      Named("forward") = fw_performances, 
      Named("backward") = bw_performances
    ),
    Named("model")=List::create( 
      Named("forward") = fw_models, 
      Named("backward") = bw_models 
    )
  );

  return(res);

}

#endif
 
