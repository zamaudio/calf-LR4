/* Calf DSP Library
 * Biquad filters
 * Copyright (C) 2001-2007 Krzysztof Foltman
 *
 * Most of code in this file is based on freely
 * available other work of other people (filter equations).
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301  USA
 */
#ifndef __CALF_BIQUAD_H
#define __CALF_BIQUAD_H

#include <complex>
#include "primitives.h"

#define IMAG (std::complex<double>(0,1))
#define FILTER_TYPE_HP 1
#define FILTER_TYPE_LP 2
#define MAXPZ 512 

namespace dsp {

/**
 * Coefficients for two-pole two-zero filter, for floating point values,
 * plus a bunch of functions to set them to typical values.
 * 
 * Coefficient calculation is based on famous Robert Bristow-Johnson's equations,
 * except where it's not. 
 * The coefficient calculation is NOT mine, the only exception is the lossy 
 * optimization in Zoelzer and rbj HP filter code.
 * 
 * See http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt for reference.
 * 
 * don't use this for integers because it won't work
 */
template<class Coeff = float>
class biquad_coeffs
{
public:
    // filter coefficients
    Coeff a0, a1, a2, b1, b2, gain;
    typedef std::complex<double> cfloat;

    biquad_coeffs()
    {
        set_null();
    }
    
    inline void set_null()
    {
        a0 = 1.0;
        b1 = b2 = a1 = a2 = 0.f;
    }
    
    /** Lowpass filter based on Robert Bristow-Johnson's equations
     * Perhaps every synth code that doesn't use SVF uses these
     * equations :)
     * @param fc     resonant frequency
     * @param q      resonance (gain at fc)
     * @param sr     sample rate
     * @param gain   amplification (gain at 0Hz)
     */
    inline void set_lp_rbj(float fc, float q, float sr, float gain = 1.0)
    {
        float omega=(float)(2*M_PI*fc/sr);
        float sn=sin(omega);
        float cs=cos(omega);
        float alpha=(float)(sn/(2*q));
        float inv=(float)(1.0/(1.0+alpha));

        a2 = a0 =  (float)(gain*inv*(1 - cs)*0.5f);
        a1 =  a0 + a0;
        b1 =  (float)(-2*cs*inv);
        b2 =  (float)((1 - alpha)*inv);
    }

    // different lowpass filter, based on Zoelzer's equations, modified by
    // me (kfoltman) to use polynomials to approximate tangent function
    // not very accurate, but perhaps good enough for synth work :)
    // odsr is "one divided by samplerate"
    // from how it looks, it perhaps uses bilinear transform - but who knows :)
    inline void set_lp_zoelzer(float fc, float q, float odsr, float gain=1.0)
    {
        Coeff omega=(Coeff)(M_PI*fc*odsr);
        Coeff omega2=omega*omega;
        Coeff K=omega*(1+omega2*omega2*Coeff(1.0/1.45));
        Coeff KK=K*K;
        Coeff QK=q*(KK+1.f);
        Coeff iQK=1.0f/(QK+K);
        Coeff inv=q*iQK;
        b2 =  (Coeff)(iQK*(QK-K));
        b1 =  (Coeff)(2.f*(KK-1.f)*inv);
        a2 = a0 =  (Coeff)(inv*gain*KK);
        a1 =  a0 + a0;
    }

    /** Highpass filter based on Robert Bristow-Johnson's equations
     * @param fc     resonant frequency
     * @param q      resonance (gain at fc)
     * @param sr     sample rate
     * @param gain   amplification (gain at sr/2)
     */
    inline void set_hp_rbj(float fc, float q, float esr, float gain=1.0)
    {
        Coeff omega=(float)(2*M_PI*fc/esr);
        Coeff sn=sin(omega);
        Coeff cs=cos(omega);
        Coeff alpha=(float)(sn/(2*q));

        float inv=(float)(1.0/(1.0+alpha));

        a0 =  (Coeff)(gain*inv*(1 + cs)/2);
        a1 =  -2.f * a0;
        a2 =  a0;
        b1 =  (Coeff)(-2*cs*inv);
        b2 =  (Coeff)((1 - alpha)*inv);
    }

    // this replaces sin/cos with polynomial approximation
    inline void set_hp_rbj_optimized(float fc, float q, float esr, float gain=1.0)
    {
        Coeff omega=(float)(2*M_PI*fc/esr);
        Coeff sn=omega+omega*omega*omega*(1.0/6.0)+omega*omega*omega*omega*omega*(1.0/120);
        Coeff cs=1-omega*omega*(1.0/2.0)+omega*omega*omega*omega*(1.0/24);
        Coeff alpha=(float)(sn/(2*q));

        float inv=(float)(1.0/(1.0+alpha));

        a0 =  (Coeff)(gain*inv*(1 + cs)*(1.0/2.0));
        a1 =  -2.f * a0;
        a2 =  a0;
        b1 =  (Coeff)(-2*cs*inv);
        b2 =  (Coeff)((1 - alpha)*inv);
    }
    
    /** Bandpass filter based on Robert Bristow-Johnson's equations (normalized to 1.0 at center frequency)
     * @param fc     center frequency (gain at fc = 1.0)
     * @param q      =~ fc/bandwidth (not quite, but close)  - 1/Q = 2*sinh(ln(2)/2*BW*w0/sin(w0))
     * @param sr     sample rate
     * @param gain   amplification (gain at sr/2)
     */
    inline void set_bp_rbj(double fc, double q, double esr, double gain=1.0)
    {
        float omega=(float)(2*M_PI*fc/esr);
        float sn=sin(omega);
        float cs=cos(omega);
        float alpha=(float)(sn/(2*q));

        float inv=(float)(1.0/(1.0+alpha));

        a0 =  (float)(gain*inv*alpha);
        a1 =  0.f;
        a2 =  (float)(-gain*inv*alpha);
        b1 =  (float)(-2*cs*inv);
        b2 =  (float)((1 - alpha)*inv);
    }
    
    // rbj's bandreject
    inline void set_br_rbj(double fc, double q, double esr, double gain=1.0)
    {
        float omega=(float)(2*M_PI*fc/esr);
        float sn=sin(omega);
        float cs=cos(omega);
        float alpha=(float)(sn/(2*q));

        float inv=(float)(1.0/(1.0+alpha));

        a0 =  (Coeff)(gain*inv);
        a1 =  (Coeff)(-gain*inv*2*cs);
        a2 =  (Coeff)(gain*inv);
        b1 =  (Coeff)(-2*cs*inv);
        b2 =  (Coeff)((1 - alpha)*inv);
    }
    // this is mine (and, I guess, it sucks/doesn't work)
    void set_allpass(float freq, float pole_r, float sr)
    {
        float a=prewarp(freq, sr);
        float q=pole_r;
        set_bilinear(a*a+q*q, -2.0f*a, 1, a*a+q*q, 2.0f*a, 1);
    }
    /// prewarping for bilinear transform, maps given digital frequency to analog counterpart for analog filter design
    static inline float prewarp(float freq, float sr)
    {
        if (freq>sr*0.49) freq=(float)(sr*0.49);
        return (float)(tan(M_PI*freq/sr));
    }
    /// convert analog angular frequency value to digital
    static inline float unwarp(float omega, float sr)
    {
        float T = 1.0 / sr;
        return (2 / T) * atan(omega * T / 2);
    }
    /// convert analog filter time constant to digital counterpart
    static inline float unwarpf(float t, float sr)
    {
        // this is most likely broken and works by pure accident!
        float omega = 1.0 / t;
        omega = unwarp(omega, sr);
        // I really don't know why does it have to be M_PI and not 2 * M_PI!
        float f = M_PI / omega;
        return f / sr;
    }
    /// set digital filter parameters based on given analog filter parameters
    void set_bilinear(float aa0, float aa1, float aa2, float ab0, float ab1, float ab2)
    {
        float q=(float)(1.0/(ab0+ab1+ab2));
        a0 = (aa0+aa1+aa2)*q;
        a1 = 2*(aa0-aa2)*q;
        a2 = (aa0-aa1+aa2)*q;
        b1 = 2*(ab0-ab2)*q;
        b2 = (ab0-ab1+ab2)*q;
    }
    
    /// RBJ peaking EQ
    /// @param freq   peak frequency
    /// @param q      q (correlated to freq/bandwidth, @see set_bp_rbj)
    /// @param peak   peak gain (1.0 means no peak, >1.0 means a peak, less than 1.0 is a dip)
    inline void set_peakeq_rbj(float freq, float q, float peak, float sr)
    {
        float A = sqrt(peak);
        float w0 = freq * 2 * M_PI * (1.0 / sr);
        float alpha = sin(w0) / (2 * q);
        float ib0 = 1.0 / (1 + alpha/A);
        a1 = b1 = -2*cos(w0) * ib0;
        a0 = ib0 * (1 + alpha*A);
        a2 = ib0 * (1 - alpha*A);
        b2 = ib0 * (1 - alpha/A);
    }
    
    /// RBJ low shelf EQ - amplitication of 'peak' at 0 Hz and of 1.0 (0dB) at sr/2 Hz
    /// @param freq   corner frequency (gain at freq is sqrt(peak))
    /// @param q      q (relates bandwidth and peak frequency), the higher q, the louder the resonant peak (situated below fc) is
    /// @param peak   shelf gain (1.0 means no peak, >1.0 means a peak, less than 1.0 is a dip)
    inline void set_lowshelf_rbj(float freq, float q, float peak, float sr)
    {
        float A = sqrt(peak);
        float w0 = freq * 2 * M_PI * (1.0 / sr);
        float alpha = sin(w0) / (2 * q);
        float cw0 = cos(w0);
        float tmp = 2 * sqrt(A) * alpha;
        float b0 = 0.f, ib0 = 0.f;
        
        a0 =    A*( (A+1) - (A-1)*cw0 + tmp);
        a1 =  2*A*( (A-1) - (A+1)*cw0);
        a2 =    A*( (A+1) - (A-1)*cw0 - tmp);
        b0 =        (A+1) + (A-1)*cw0 + tmp;
        b1 =   -2*( (A-1) + (A+1)*cw0);
        b2 =        (A+1) + (A-1)*cw0 - tmp;
        
        ib0 = 1.0 / b0;
        b1 *= ib0;
        b2 *= ib0;
        a0 *= ib0;
        a1 *= ib0;
        a2 *= ib0;
    }
    
    /// RBJ high shelf EQ - amplitication of 0dB at 0 Hz and of peak at sr/2 Hz
    /// @param freq   corner frequency (gain at freq is sqrt(peak))
    /// @param q      q (relates bandwidth and peak frequency), the higher q, the louder the resonant peak (situated above fc) is
    /// @param peak   shelf gain (1.0 means no peak, >1.0 means a peak, less than 1.0 is a dip)
    inline void set_highshelf_rbj(float freq, float q, float peak, float sr)
    {
        float A = sqrt(peak);
        float w0 = freq * 2 * M_PI * (1.0 / sr);
        float alpha = sin(w0) / (2 * q);
        float cw0 = cos(w0);
        float tmp = 2 * sqrt(A) * alpha;
        float b0 = 0.f, ib0 = 0.f;
        
        a0 =    A*( (A+1) + (A-1)*cw0 + tmp);
        a1 = -2*A*( (A-1) + (A+1)*cw0);
        a2 =    A*( (A+1) + (A-1)*cw0 - tmp);
        b0 =        (A+1) - (A-1)*cw0 + tmp;
        b1 =    2*( (A-1) - (A+1)*cw0);
        b2 =        (A+1) - (A-1)*cw0 - tmp;
        
        ib0 = 1.0 / b0;
        b1 *= ib0;
        b2 *= ib0;
        a0 *= ib0;
        a1 *= ib0;
        a2 *= ib0;
    }
    
/// Butterworth digital filter coefficient calculator code
/// thanks to late Tony Fisher of University of York, Computer Science Dept
void choosepole(std::complex<double> z)
{
	if (real(z) < 0.0)
        {
		if (polemask & 1) splane.poles[splane.numpoles++] = z;
		polemask >>= 1;
	}
}

void prewarp() {
	warped_alpha1 = tan(M_PI * raw_alpha1) / M_PI;
	warped_alpha2 = tan(M_PI * raw_alpha2) / M_PI;
}

void normalize(int type) {
	double w1 = 2 * M_PI * warped_alpha1;
	if(type==FILTER_TYPE_LP) {
		for (int i = 0; i < splane.numpoles; i++) {
			splane.poles[i] = splane.poles[i] * w1;
		}
		splane.numzeros = 0;
	}

	if(type==FILTER_TYPE_HP) {
		for (int i=0; i < splane.numpoles; i++) splane.poles[i] = w1 / splane.poles[i];
		for (int i=0; i < splane.numpoles; i++) splane.zeros[i] = 0.0;	 /* also N zeros at (0,0) */
		splane.numzeros = splane.numpoles;
	}
}

void multin(std::complex<double> w, int npz, std::complex<double> coeffs[]){
	/* multiply factor (z-w) into coeffs */
	std::complex<double> nw = -w;
	for (int i = npz; i >= 1; i--) coeffs[i] = (nw * coeffs[i]) + coeffs[i-1];
	coeffs[0] = nw * coeffs[0];
}

std::complex<double> eval(std::complex<double> coeffs[], int npz, std::complex<double> z) const {
	 std::complex<double> sum = 0.0;
	 for (int i = npz; i >= 0; i--) sum = (sum * z) + coeffs[i];
	 return sum;
}

std::complex<double> evaluate(std::complex<double> topco[], int nz, std::complex<double> botco[], int np, std::complex<double> z) const
{ /* evaluate response, substituting for z */
	return eval(topco, nz, z) / eval(botco, np, z);
}

void expand(std::complex<double> pz[], int npz, std::complex<double> coeffs[]){
	int i;
	coeffs[0] = 1.0;
	for (i=0; i < npz; i++) coeffs[i+1] = 0.0;
	for (i=0; i < npz; i++) multin(pz[i], npz, coeffs);
}

void expandpoly(int type, double *xcoeffs, double *ycoeffs, double *final_gain) {
	std::complex<double> topcoeffs[3], botcoeffs[3]; int i;
	expand(zplane.zeros, zplane.numzeros, topcoeffs);
	expand(zplane.poles, zplane.numpoles, botcoeffs);
	dc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, 1.0);
	double theta = 2.0 * M_PI * raw_alpha1;
	fc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, exp(IMAG*theta));
	hf_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, -1.0);
	for (i = 0; i <= zplane.numzeros; i++) xcoeffs[i] = real(topcoeffs[i]) / real(botcoeffs[zplane.numpoles]);
	for (i = 0; i <= zplane.numpoles; i++) ycoeffs[i] = -real(botcoeffs[i]) / real(botcoeffs[zplane.numpoles]);

	if(type==FILTER_TYPE_HP) {
		*final_gain=real(hf_gain);
	}
	if(type==FILTER_TYPE_LP) {
		*final_gain=real(dc_gain);
	}
}

void compute_s() {
	for (int i = 0; i < 2*order; i++)
	{
		double theta = (order & 1) ? (i*M_PI) / order : ((i+0.5)*M_PI) / order;
		choosepole(exp(IMAG*theta));
	}
}

void compute_z_blt() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by bilinear transform */
 {
	int i;
	zplane.numpoles = splane.numpoles;
	zplane.numzeros = splane.numzeros;
	for (i=0; i < zplane.numpoles; i++) zplane.poles[i] = (2.0+splane.poles[i])/(2.0-splane.poles[i]);
	for (i=0; i < zplane.numzeros; i++) zplane.zeros[i] = (2.0+splane.zeros[i])/(2.0-splane.zeros[i]);
	while (zplane.numzeros < zplane.numpoles) zplane.zeros[zplane.numzeros++] = -1.0;
}

int order;
unsigned int polemask;
double xcoeffs[4], ycoeffs[4];
double raw_alpha1, raw_alpha2;
cfloat dc_gain, fc_gain, hf_gain;
double final_gain;
double warped_alpha1, warped_alpha2;

struct pzrep
{ 
	std::complex<double> poles[3], zeros[3];
        int numpoles, numzeros;
};
pzrep splane, zplane;
      
    inline void set_hp_lr4(float freq, float sr)
    {
    	int type = 1;
	raw_alpha1 = freq / sr;
	raw_alpha2 = raw_alpha1;
        polemask = ~0;
	order = 2;
	splane.numpoles = 0;
	compute_s();
	prewarp();
	normalize(type);
	compute_z_blt();
	expandpoly(type, &xcoeffs[0], &ycoeffs[0], &final_gain);
	a0 = xcoeffs[0];
	a1 = xcoeffs[1];
	a2 = xcoeffs[2];
	b1 = ycoeffs[0];
	b2 = ycoeffs[1];
	gain = final_gain;
	printf("hp: a0=%f a1=%f a2=%f b1=%f b2=%f gain=%f\n",a0,a1,a2,b1,b2,gain);
    }

    inline void set_lp_lr4(float freq, float sr)
    {
    	int type = 2;
	raw_alpha1 = freq / sr;
	raw_alpha2 = raw_alpha1;
        polemask = ~0;
	order = 2;
	splane.numpoles = 0;
	compute_s();
	prewarp();
	normalize(type);
	compute_z_blt();
	expandpoly(type, &xcoeffs[0], &ycoeffs[0], &final_gain);
	a0 = xcoeffs[0];
	a1 = xcoeffs[1];
	a2 = xcoeffs[2];
	b1 = ycoeffs[0];
	b2 = ycoeffs[1];
	gain = final_gain;
	printf("lp: a0=%f a1=%f a2=%f b1=%f b2=%f gain=%f\n",a0,a1,a2,b1,b2,gain);
    }

    /// copy coefficients from another biquad
    template<class U>
    inline void copy_coeffs(const biquad_coeffs<U> &src)
    {
        a0 = src.a0;
        a1 = src.a1;
        a2 = src.a2;
        b1 = src.b1;
        b2 = src.b2;
    }
    

float freq_gain(float freq, float sr) const 
{ 

    typedef std::complex<double> cfloat;
    int nzeros = 2;
    int npoles = 2;
    cfloat topcos[MAXPZ+1], botcos[MAXPZ+1]; 
    topcos[0] = cfloat(a0, 0.0);
    topcos[1] = cfloat(a1, 0.0);
    topcos[2] = cfloat(a2, 0.0);	
    botcos[0] = cfloat(b1, 0.0);
    botcos[1] = cfloat(b2, 0.0);
    
    double theta = 2.0 * M_PI * freq / sr;

    cfloat z = exp(cfloat(0.0, theta));
    cfloat fr = evaluate(topcos, nzeros, botcos, npoles, z);
/*
    double theta = 2.0 * M_PI * freq / sr;
    cfloat zinv  = exp(cfloat(0.0, -theta));
    cfloat fr = cfloat(1.0, 0.0) / ((cfloat(a0, 0.0) + cfloat(a1, 0.0) * zinv + cfloat(a2, 0.0) * zinv * zinv) / (cfloat(1.0, 0.0) + cfloat(b1, 0.0) * zinv + cfloat(b2, 0.0) * zinv * zinv));
*/
    float mag = std::abs(fr);
    //printf("(%.2f,%.2f) ", freq, mag);
    return mag;
}

    /// Return the filter's gain at frequency freq
    /// @param freq   Frequency to look up
    /// @param sr     Filter sample rate (used to convert frequency to angular frequency)
    float freq_gain_old(float freq, float sr) const
    {
        typedef std::complex<double> cfloat;
        freq *= 2.0 * M_PI / sr;
        cfloat z = 1.0 / exp(cfloat(0.0, freq));
        
        return std::abs(h_z(z));
    }

    /// Return H(z) the filter's gain at frequency freq
    /// @param z   Z variable (e^jw)
    cfloat h_z(const cfloat &z) const
    {
        
        return (cfloat(a0) + double(a1) * z + double(a2) * z*z) / (cfloat(1.0) + double(b1) * z + double(b2) * z*z);
    }
    
};

/**
 * Two-pole two-zero filter, for floating point values.
 * Uses "traditional" Direct I form (separate FIR and IIR halves).
 * don't use this for integers because it won't work
 */
template<class Coeff = float, class T = float>
struct biquad_d1: public biquad_coeffs<Coeff>
{
    using biquad_coeffs<Coeff>::a0;
    using biquad_coeffs<Coeff>::a1;
    using biquad_coeffs<Coeff>::a2;
    using biquad_coeffs<Coeff>::b1;
    using biquad_coeffs<Coeff>::b2;
    /// input[n-1]
    T x1; 
    /// input[n-2]
    T x2; 
    /// output[n-1]
    T y1; 
    /// output[n-2]
    T y2; 
    /// Constructor (initializes state to all zeros)
    biquad_d1()
    {
        reset();
    }
    /// direct I form with four state variables
    inline T process(T in)
    {
        T out = in * a0 + x1 * a1 + x2 * a2 - y1 * b1 - y2 * b2;
        x2 = x1;
        y2 = y1;
        x1 = in;
        y1 = out;
        return out;
    }
    
    /// direct I form with zero input
    inline T process_zeroin()
    {
        T out = - y1 * b1 - y2 * b2;
        y2 = y1;
        y1 = out;
        return out;
    }
    
    /// simplified version for lowpass case with two zeros at -1
    inline T process_lp(T in)
    {
        T out = a0*(in + x1 + x1 + x2) - y1 * b1 - y2 * b2;
        x2 = x1;
        y2 = y1;
        x1 = in;
        y1 = out;
        return out;
    }
    /// Sanitize (set to 0 if potentially denormal) filter state
    inline void sanitize() 
    {
        dsp::sanitize(x1);
        dsp::sanitize(y1);
        dsp::sanitize(x2);
        dsp::sanitize(y2);
    }
    /// Reset state variables
    inline void reset()
    {
        dsp::zero(x1);
        dsp::zero(y1);
        dsp::zero(x2);
        dsp::zero(y2);
    }
    inline bool empty() const {
        return (y1 == 0.f && y2 == 0.f);
    }
    
};
    
/**
 * Two-pole two-zero filter, for floating point values.
 * Uses slightly faster Direct II form (combined FIR and IIR halves).
 * However, when used with wildly varying coefficients, it may 
 * make more zipper noise than Direct I form, so it's better to
 * use it when filter coefficients are not changed mid-stream.
 */
template<class Coeff = float, class T = float>
struct biquad_d2: public biquad_coeffs<Coeff>
{
    using biquad_coeffs<Coeff>::a0;
    using biquad_coeffs<Coeff>::a1;
    using biquad_coeffs<Coeff>::a2;
    using biquad_coeffs<Coeff>::b1;
    using biquad_coeffs<Coeff>::b2;
    /// state[n-1]
    float w1; 
    /// state[n-2]
    float w2; 
    /// Constructor (initializes state to all zeros)
    biquad_d2()
    {
        reset();
    }
    /// direct II form with two state variables
    inline T process(T in)
    {
        dsp::sanitize_denormal(in);
        dsp::sanitize(in);
        dsp::sanitize(w1);
        dsp::sanitize(w2);

        T tmp = in - w1 * b1 - w2 * b2;
        T out = tmp * a0 + w1 * a1 + w2 * a2;
        w2 = w1;
        w1 = tmp;
        return out;
    }
    
    // direct II form with two state variables, lowpass version
    // interesting fact: this is actually slower than the general version!
    inline T process_lp(T in)
    {
        T tmp = in - w1 * b1 - w2 * b2;
        T out = (tmp  + w2 + w1* 2) * a0;
        w2 = w1;
        w1 = tmp;
        return out;
    }

    /// Is the filter state completely silent? (i.e. set to 0 by sanitize function)
    inline bool empty() const {
        return (w1 == 0.f && w2 == 0.f);
    }
    
    
    /// Sanitize (set to 0 if potentially denormal) filter state
    inline void sanitize() 
    {
        dsp::sanitize(w1);
        dsp::sanitize(w2);
    }
    
    /// Reset state variables
    inline void reset()
    {
        dsp::zero(w1);
        dsp::zero(w2);
    }
};

/**
 * Two-pole two-zero filter, for floating point values.
 * Uses "traditional" Direct I form (separate FIR and IIR halves).
 * don't use this for integers because it won't work
 */
template<class Coeff = float, class T = float>
struct biquad_d1_lerp: public biquad_coeffs<Coeff>
{
    using biquad_coeffs<Coeff>::a0;
    using biquad_coeffs<Coeff>::a1;
    using biquad_coeffs<Coeff>::a2;
    using biquad_coeffs<Coeff>::b1;
    using biquad_coeffs<Coeff>::b2;
    Coeff a0cur, a1cur, a2cur, b1cur, b2cur;
    Coeff a0delta, a1delta, a2delta, b1delta, b2delta;
    /// input[n-1]
    T x1; 
    /// input[n-2]
    T x2; 
    /// output[n-1]
    T y1; 
    /// output[n-2]
    T y2; 
    /// Constructor (initializes state to all zeros)
    biquad_d1_lerp()
    {
        reset();
    }
    #define _DO_COEFF(coeff) coeff##delta = (coeff - coeff##cur) * (frac)
    void big_step(Coeff frac)
    {
        _DO_COEFF(a0);
        _DO_COEFF(a1);
        _DO_COEFF(a2);
        _DO_COEFF(b1);
        _DO_COEFF(b2);
    }
    #undef _DO_COEFF
    /// direct I form with four state variables
    inline T process(T in)
    {
        T out = in * a0cur + x1 * a1cur + x2 * a2cur - y1 * b1cur - y2 * b2cur;
        x2 = x1;
        y2 = y1;
        x1 = in;
        y1 = out;
        a0cur += a0delta;
        a1cur += a1delta;
        a2cur += a2delta;
        b1cur += b1delta;
        b2cur += b2delta;
        return out;
    }
    
    /// direct I form with zero input
    inline T process_zeroin()
    {
        T out = - y1 * b1 - y2 * b2;
        y2 = y1;
        y1 = out;
        b1cur += b1delta;
        b2cur += b2delta;
        return out;
    }
    
    /// simplified version for lowpass case with two zeros at -1
    inline T process_lp(T in)
    {
        T out = a0*(in + x1 + x1 + x2) - y1 * b1 - y2 * b2;
        x2 = x1;
        y2 = y1;
        x1 = in;
        y1 = out;
        return out;
    }
    /// Sanitize (set to 0 if potentially denormal) filter state
    inline void sanitize() 
    {
        dsp::sanitize(x1);
        dsp::sanitize(y1);
        dsp::sanitize(x2);
        dsp::sanitize(y2);
        dsp::sanitize(a0cur);
        dsp::sanitize(a1cur);
        dsp::sanitize(a2cur);
        dsp::sanitize(b1cur);
        dsp::sanitize(b2cur);
    }
    /// Reset state variables
    inline void reset()
    {
        dsp::zero(x1);
        dsp::zero(y1);
        dsp::zero(x2);
        dsp::zero(y2);
        dsp::zero(a0cur);
        dsp::zero(a1cur);
        dsp::zero(a2cur);
        dsp::zero(b1cur);
        dsp::zero(b2cur);
    }
    inline bool empty() {
        return (y1 == 0.f && y2 == 0.f);
    }
    
};
    

/**
 * Damien's LR4 filter routine:
 * Two poles, two IIR Butterworth filters cascaded Q=1/sqrt(2).
 */
template<class Coeff = float, class T = float>
struct biquad_lr4: public biquad_coeffs<Coeff>
{
    using biquad_coeffs<Coeff>::a0;
    using biquad_coeffs<Coeff>::a1;
    using biquad_coeffs<Coeff>::a2;
    using biquad_coeffs<Coeff>::b1;
    using biquad_coeffs<Coeff>::b2;
    using biquad_coeffs<Coeff>::gain;
    /// Filter 1
    /// input
    T x0; 
    /// input
    T x1; 
    /// input
    T x2;
    /// output
    T y0; 
    /// output
    T y1; 
    /// output
    T y2;

    /// Filter 2 (cascaded)
    /// input
    T xx0; 
    /// input
    T xx1; 
    /// input
    T xx2;
    /// output
    T yy0; 
    /// output
    T yy1; 
    /// output
    T yy2;
    
    /// Constructor (initializes state to all zeros)
    biquad_lr4()
    {
        reset();
    }
    /// process with twelve state variables
    inline T process(T in)
    {
        x0 = x1;
        x1 = x2;
        x2 = in / gain;
        y0 = y1;
        y1 = y2;
        y2 = a0 * x0 + a1 * x1 + a2 * x2 + b1 * y0 + b2 * y1;
        
	T temp = y2;

	xx0 = xx1;
        xx1 = xx2;
        xx2 = temp / gain;
        yy0 = yy1;
        yy1 = yy2;
        yy2 = a0 * xx0 + a1 * xx1 + a2 * xx2 + b1 * yy0 + b2 * yy1;
	
	T out = yy2;
        return out;
    }
    
    /// zero input
    inline T process_zeroin()
    {
        T zero = 0.f;
	return process(zero);
    }
    
    /// simplified version for lowpass case with two zeros at -1
    inline T process_lp(T in)
    {
        return process(in);
    }

    /// Sanitize (set to 0 if potentially denormal) filter state
    inline void sanitize() 
    {
/*
	dsp::sanitize(x0);
        dsp::sanitize(y0);
        dsp::sanitize(x1);
        dsp::sanitize(y1);
        dsp::sanitize(x2);
        dsp::sanitize(y2);
*/
    }
    /// Reset state variables
    inline void reset()
    {
    	dsp::zero(x0);
	dsp::zero(y0);
        dsp::zero(x1);
        dsp::zero(y1);
        //dsp::zero(x2);
        //dsp::zero(y2);
    }
    inline bool empty() const {
        return (y0 == 0.f && y1 == 0.f && y2 == 0.f);
    }
    
};



/// Compose two filters in series
template<class F1, class F2>
class filter_compose {
public:
    typedef std::complex<float> cfloat;
    F1 f1;
    F2 f2;
public:
    float process(float value) {
        return f2.process(f1.process(value));
    }
    
    inline cfloat h_z(const cfloat &z) const {
        return f1.h_z(z) * f2.h_z(z);
    }
    
    /// Return the filter's gain at frequency freq
    /// @param freq   Frequency to look up
    /// @param sr     Filter sample rate (used to convert frequency to angular frequency)
    float freq_gain(float freq, float sr) const
    {
        typedef std::complex<double> cfloat;
        freq *= 2.0 * M_PI / sr;
        cfloat z = 1.0 / exp(cfloat(0.0, freq));
        
        return std::abs(h_z(z));
    }
    
    void sanitize() {
        f1.sanitize();
        f2.sanitize();
    }
};

/// Compose two filters in parallel
template<class F1, class F2>
class filter_sum {
public:
    typedef std::complex<double> cfloat;
    F1 f1;
    F2 f2;
public:
    float process(float value) {
        return f2.process(value) + f1.process(value);
    }
    
    inline cfloat h_z(const cfloat &z) const {
        return f1.h_z(z) + f2.h_z(z);
    }
    
    /// Return the filter's gain at frequency freq
    /// @param freq   Frequency to look up
    /// @param sr     Filter sample rate (used to convert frequency to angular frequency)
    float freq_gain(float freq, float sr) const
    {
        typedef std::complex<double> cfloat;
        freq *= 2.0 * M_PI / sr;
        cfloat z = 1.0 / exp(cfloat(0.0, freq));
        
        return std::abs(h_z(z));
    }
    
    void sanitize() {
        f1.sanitize();
        f2.sanitize();
    }
};

};

#endif
