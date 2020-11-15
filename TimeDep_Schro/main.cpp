//
//  main.cpp
//  TimeDep_Schro
//
//  Created by Ryan Hill on 11/10/20.
//
//  Program to calculate the time dependent propogation of an
//  electron wave packet through a potential barrier.

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include "arrayt.hpp"

typedef std::complex<double> CMPLX;
typedef arrayt<CMPLX> arrayc;

// Parameters
const int L = 500; // length of domain in space (Angstroms)
const int Tf = 500; // length of domain in time (1e-12 seconds)
const int s = 10; // Gaussian wavefunction width (Angstroms)
const double k0 = 1; // Gaussian wavefunction average wavenumber (inverse Angstroms)
const double V1 = 2.0; // potential energy coefficient (eV)
const double x0 = 0.5*L; // electron wavepacket initial position (Angstroms)
const int wx = 5; // width of potential peaks (Angstroms)
const int Nx = 500; // number of grid points (dimensionless)
// const int Nt = 500; // number of time steps (dimensionless)
const double dx = L/Nx; // sampling size in space (Angstrom)
const double dt = 1; // sampling size in time (1e-12seconds)
const double w = 2.0*dx*dx/dt; // paramter found in FD form of Schro eq. (Angstroms^2/second)
const double hbar = 6.5821e-16; // Planck's constant (eV-sec)
const CMPLX I(0.0, 1.0); // imaginary unit
const double c0 = 3.801; // hbar^2/2m for electron (eV-A^2)
const CMPLX c1 = I*w*hbar/c0; // 2mwi/hbar
const CMPLX c2 = dx*dx/c0; // 2mdx^2/hbar^2

std::ofstream outfile_psi_re0 ("/Users/ryanhill/Desktop/TimeDep_Schro/hw07-psi-re0.dat");
std::ofstream outfile_psi_im0 ("/Users/ryanhill/Desktop/TimeDep_Schro/hw07-psi-im0.dat");
std::ofstream outfile_psi_sq0 ("/Users/ryanhill/Desktop/TimeDep_Schro/hw07-psi-sq0.dat");
std::ofstream outfile_psi_sq ("/Users/ryanhill/Desktop/TimeDep_Schro/hw07-psi-sq25.dat");
std::ofstream outfile_v ("/Users/ryanhill/Desktop/TimeDep_Schro/hw07-v.dat");
std::ofstream outfile_x ("/Users/ryanhill/Desktop/TimeDep_Schro/hw07-x.dat");

// Initial complex valued Gaussian wavefunction for an electron.
CMPLX calc_psi0(double x) {
    CMPLX term1 = (x-0.3*L)/(s);
    CMPLX term2 = I*x*k0;
    return exp(-(term1*term1)+term2);
}

// Potential energy modeling a one dimensional crystal surface with periodic peaks
double calc_V(double x) {
    if (x>x0) { return V1*(0.75-cos((x-x0)/wx)); }
    else { return 0; }
}

// Implements the implicit Crank-Nicolson method
void cn_method(arrayc &psi, arrayt<double> &V) {
    
    // define four 1D arrays
    arrayc a(Nx);
    arrayc b(Nx);
    arrayc c(Nx);
    arrayc d(Nx);
    
    for (int j=0; j<Nx; j++) {
        a(j) = 1;
        b(j) = c1-2.0-c2*V(j+1);
        c(j) = 1;
        d(j) = -psi(j)+(c1+2.0+c2*V(j+1))*psi(j+1)-psi(j+2);
    }
    
    // exclude first member of left diagonal and last member of right diagonal
    a(0) = 0.0;
    c(Nx-1) = 0.0;
    
    // update coefficients in first row
    c(0) = c(0)/b(0);
    d(0) = d(0)/b(0);
    
    d(Nx-1) = 0.0; // boundary condition
    
    // setup tri-diagonal matrix by reducing to upper diagonal form and then using backsubstitution
    for (int n=0; n<Tf; n++) {
        
        for (int j=1; j<Nx-2; j++) {
            c(j) = c(j)/(b(j)-a(j)*c(j-1));
        }
        
        for (int j=1; j<Nx-1; j++) {
            d(j) = (d(j)-a(j)*d(j-1))/(b(j)-a(j)*c(j-1));
        }
        
        psi(Nx-1)=d(Nx-1);
        
        for (int j=Nx-2; j>=0; j--) {
            psi(j) = d(j)-c(j)*psi(j+1);
        }
    }
}

int main() {
    
    // define arrays for wavefunction, potential
    arrayc psi(Nx);
    arrayt<double> V(Nx);
    
    // initialize psi and V, populate x values
    for (int j=0; j<Nx; j++) {
        psi(j) = calc_psi0(j*dx);
        V(j) = calc_V(j*dx);
        outfile_psi_re0 << psi(j).real() << std::endl;
        outfile_psi_im0 << psi(j).imag() << std::endl;
        outfile_psi_sq0 << abs((psi(j)*psi(j)).real()) << std::endl;
        outfile_v << V(j) << std::endl;
    }
    
    cn_method(psi, V);
    
    // store final values in outfiles
    for (int j=0; j<Nx; j++) {
        outfile_psi_sq << abs((psi(j)*psi(j)).real()) << std::endl;
        outfile_x << j*dx << std::endl;
    }
    outfile_psi_re0 << std::endl; // empty line
    outfile_psi_im0 << std::endl;
    outfile_psi_sq0 << std::endl;
    outfile_psi_sq << std::endl;
    outfile_v << std::endl;
    outfile_x << std::endl;
        
    outfile_psi_re0.close();
    outfile_psi_im0.close();
    outfile_psi_sq0.close();
    outfile_psi_sq.close();
    outfile_v.close();
    outfile_x.close();
    
    std::cout << "finished for t=" << Tf/200.0 << "e-14 sec" << endl;
    return 0;
    
} // end main
