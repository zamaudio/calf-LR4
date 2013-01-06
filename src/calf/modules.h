/* Calf DSP plugin pack
 * Assorted plugins
 *
 * Copyright (C) 2001-2010 Krzysztof Foltman, Markus Schmidt, Thor Harald Johansen and others
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
 * Boston, MA 02111-1307, USA.
 */
#ifndef CALF_MODULES_H
#define CALF_MODULES_H

#include <assert.h>
#include <fftw3.h>
#include <limits.h>
#include "biquad.h"
#include "inertia.h"
#include "audio_fx.h"
#include "giface.h"
#include "metadata.h"
#include "loudness.h"

#include <iostream>
#include <cmath>
#include <complex>

namespace calf_plugins {
#define IMAG (std::complex<double>(0,1))

#define FILTER_TYPE_HP 1

#define FILTER_TYPE_LP 2

#define LP1L 0
#define HP1L 1
#define LP2L 2
#define HP2L 3
#define LP1R 4
#define HP1R 5
#define LP2R 6
#define HP2R 7


struct pzrep

{ std::complex<double> poles[3], zeros[3];
    int numpoles, numzeros;
 };

// Butterworth digital filter coefficient calculator code
// thanks to late Tony Fisher of University of York, Computer Science Dept

class ButterworthCalculator {

 pzrep splane, zplane;
 int order;


 int type; // HP or LP
 double raw_alpha1, raw_alpha2;
 std::complex<double> dc_gain, fc_gain, hf_gain;
 double final_gain;
 double warped_alpha1, warped_alpha2;
 unsigned int polemask;
 double xcoeffs[4], ycoeffs[4];

void choosepole(std::complex<double> z)
{
	if (real(z) < 0.0)
        {
		if (polemask & 1) splane.poles[splane.numpoles++] = z;
		polemask >>= 1;
	}
}

void prewarp(){
	warped_alpha1 = tan(M_PI * raw_alpha1) / M_PI;
	warped_alpha2 = tan(M_PI * raw_alpha2) / M_PI;
}

void normalize()	{
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



std::complex<double> eval(std::complex<double> coeffs[], int npz, std::complex<double> z) {
	 std::complex<double> sum = 0.0;
	 for (int i = npz; i >= 0; i--) sum = (sum * z) + coeffs[i];
	 return sum;
}



std::complex<double> evaluate(std::complex<double> topco[], int nz, std::complex<double> botco[], int np, std::complex<double> z)
{ /* evaluate response, substituting for z */
	return eval(topco, nz, z) / eval(botco, np, z);
}


void expand(std::complex<double> pz[], int npz, std::complex<double> coeffs[]){
	int i;
	coeffs[0] = 1.0;
	for (i=0; i < npz; i++) coeffs[i+1] = 0.0;
	for (i=0; i < npz; i++) multin(pz[i], npz, coeffs);
}

void expandpoly() {
	std::complex<double> topcoeffs[3], botcoeffs[3]; int i;
	expand(zplane.zeros, zplane.numzeros, topcoeffs);
	expand(zplane.poles, zplane.numpoles, botcoeffs);
	dc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, 1.0);
	double theta = 2* M_PI * raw_alpha1;
	fc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, exp(IMAG*theta));
	hf_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, -1.0);
	for (i = 0; i <= zplane.numzeros; i++) xcoeffs[i] = real(topcoeffs[i]) / real(botcoeffs[zplane.numpoles]);
	for (i = 0; i <= zplane.numpoles; i++) ycoeffs[i] = -real(botcoeffs[i]) / real(botcoeffs[zplane.numpoles]);

	if(type==FILTER_TYPE_HP) {
		final_gain=real(hf_gain);
	}
	if(type==FILTER_TYPE_LP) {
		final_gain=real(dc_gain);
	}
}



void compute_s() {
	for (int i = 0; i < 2*order; i++)
	{
		double theta = (order & 1) ? (i*M_PI) / order : ((i+0.5)*M_PI) / order;
		choosepole(exp(IMAG*theta));
	}
}

void printrecurrence(){
	 std::cout << xcoeffs[0]<<" "<< xcoeffs[1]<<" " << xcoeffs[2]<<std::endl;
	 std::cout << ycoeffs[0]<<" "<< ycoeffs[1]<<std::endl;
	 std::cout << final_gain<<std::endl;
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

public :
void run(int ftype, double edge, double srate, double & a1, double & a2, double & a3, double & b1, double & b2, double & gain) {
	type=ftype;
	polemask = ~0;
	order=2;
	splane.numpoles = 0;
	raw_alpha1=edge/srate;
	raw_alpha2 = raw_alpha1;
	compute_s();
	prewarp();
	normalize();
	compute_z_blt();
	expandpoly();

	a1 = xcoeffs[0];

	a2 = xcoeffs[1];

	a3 = xcoeffs[2];

	b1 = ycoeffs[0];

	b2 = ycoeffs[1];

	gain = final_gain;

//	printrecurrence();
}

};



// linkwitz-riley active crossover
class xover_audio_module: public audio_module<xover_metadata>
{
    using audio_module<xover_metadata>::ins;
    using audio_module<xover_metadata>::outs;
    using audio_module<xover_metadata>::params;

private:
    float xv[4][3], yv[4][3], xv2[4][3], yv2[4][3];
    double a1[8],a2[8],a3[8],b1[8],b2[8],gain[8];
    float chain;
    ButterworthCalculator b ;
public:
    uint32_t srate;
    //static parameter_properties param_props[];

    void params_changed() {
       b.run(FILTER_TYPE_LP, *params[par_xover1],48000.0,a1[LP1L],a2[LP1L],a3[LP1L],b1[LP1L],b2[LP1L],gain[LP1L]);
	b.run(FILTER_TYPE_HP, *params[par_xover1],48000.0,a1[HP1L],a2[HP1L],a3[HP1L],b1[HP1L],b2[HP1L],gain[HP1L]);
	b.run(FILTER_TYPE_LP, *params[par_xover2],48000.0,a1[LP2L],a2[LP2L],a3[LP2L],b1[LP2L],b2[LP2L],gain[LP2L]);
	b.run(FILTER_TYPE_HP, *params[par_xover2],48000.0,a1[HP2L],a2[HP2L],a3[HP2L],b1[HP2L],b2[HP2L],gain[HP2L]);
    }

    xover_audio_module() {
	for(int i=0;i<4;++i) {
		for(int j=0;j<3;++j) {
			xv[i][j]=0.0;
			yv[i][j]=0.0;
			xv2[i][j]=0.0;
			yv2[i][j]=0.0;
		}
	}
        //*params[par_xover1]=400.0;
	//*params[par_xover2]=1400.0;
	//get coefficients for xover network at 400Hz & 1400Hz, samplerate 48kHz
	for (int  i = 0; i<3; i+=2) {

	    b.run(FILTER_TYPE_LP, 400.0+500.0*i,48000.0,a1[i],a2[i],a3[i],b1[i],b2[i],gain[i]);
	    b.run(FILTER_TYPE_HP, 400.0+500.0*i,48000.0,a1[i+1],a2[i+1],a3[i+1],b1[i+1],b2[i+1],gain[i+1]);
	}
    }

    uint32_t process(uint32_t offset, uint32_t numsamples, uint32_t inputs_mask, uint32_t outputs_mask) {
        params_changed();

	if (!inputs_mask) return 0;

	numsamples += offset;
        for (uint32_t i = offset; i < numsamples; ++i) {

	     outs[0][i] = runfilter(0,ins[0][i],a1[LP1L],a2[LP1L],a3[LP1L],b1[LP1L],b2[LP1L],gain[LP1L]);
            chain      = runfilter(1,ins[0][i],a1[HP1L],a2[HP1L],a3[HP1L],b1[HP1L],b2[HP1L],gain[HP1L]);
            outs[1][i] = runfilter(2,chain    ,a1[LP2L],a2[LP2L],a3[LP2L],b1[LP2L],b2[LP2L],gain[LP2L]);
	     outs[2][i] = runfilter(3,ins[0][i],a1[HP2L],a2[HP2L],a3[HP2L],b1[HP2L],b2[HP2L],gain[HP2L]);

	}

        return inputs_mask;
    }

    float runfilter(int n,float input,float a1,float a2,float a3,float b1,float b2,float gain) {
        float intermediate = 0;

        /* LR-4 Digital Filter */
        xv[n][0] = xv[n][1]; xv[n][1] = xv[n][2];
        xv[n][2] = input / gain;
        yv[n][0] = yv[n][1]; yv[n][1] = yv[n][2];
        yv[n][2] =   (a1*xv[n][0] + a3*xv[n][2]) + a2*xv[n][1]
                     + (b1*yv[n][0]) + (b2*yv[n][1]);
        intermediate = yv[n][2];

        xv2[n][0] = xv2[n][1]; xv2[n][1] = xv2[n][2];
        xv2[n][2] = intermediate / gain;
        yv2[n][0] = yv2[n][1]; yv2[n][1] = yv2[n][2];
        yv2[n][2] =   (a1*xv2[n][0] + a3*xv2[n][2]) + a2*xv2[n][1]
                     + (b1*yv2[n][0]) + (b2*yv2[n][1]);
        return yv2[n][2];
    }

};


struct ladspa_plugin_info;

class reverb_audio_module: public audio_module<reverb_metadata>
{
public:    
    dsp::reverb reverb;
    dsp::simple_delay<16384, dsp::stereo_sample<float> > pre_delay;
    dsp::onepole<float> left_lo, right_lo, left_hi, right_hi;
    uint32_t srate;
    dsp::gain_smoothing amount, dryamount;
    int predelay_amt;
    float meter_wet, meter_out;
    uint32_t clip;
    
    void params_changed();
    uint32_t process(uint32_t offset, uint32_t numsamples, uint32_t inputs_mask, uint32_t outputs_mask);
    void activate();
    void set_sample_rate(uint32_t sr);
    void deactivate();
};

class vintage_delay_audio_module: public audio_module<vintage_delay_metadata>
{
public:    
    // 1MB of delay memory per channel... uh, RAM is cheap
    enum { MAX_DELAY = 262144, ADDR_MASK = MAX_DELAY - 1 };
    enum { MIXMODE_STEREO, MIXMODE_PINGPONG, MIXMODE_LR, MIXMODE_RL }; 
    float buffers[2][MAX_DELAY];
    int bufptr, deltime_l, deltime_r, mixmode, medium, old_medium;
    /// number of table entries written (value is only important when it is less than MAX_DELAY, which means that the buffer hasn't been totally filled yet)
    int age;
    
    dsp::gain_smoothing amt_left, amt_right, fb_left, fb_right, dry, chmix;
    
    dsp::biquad_d2<float> biquad_left[2], biquad_right[2];
    
    uint32_t srate;
    
    vintage_delay_audio_module();
    
    void params_changed();
    void activate();
    void deactivate();
    void set_sample_rate(uint32_t sr);
    void calc_filters();
    uint32_t process(uint32_t offset, uint32_t numsamples, uint32_t inputs_mask, uint32_t outputs_mask);
    
    long _tap_avg;
    long _tap_last;
};

template<typename FilterClass, typename Metadata>
class filter_module_with_inertia: public audio_module<Metadata>, public FilterClass
{
public:
    /// These are pointers to the ins, outs, params arrays in the main class
    typedef filter_module_with_inertia inertia_filter_module;
    using audio_module<Metadata>::ins;
    using audio_module<Metadata>::outs;
    using audio_module<Metadata>::params;
    
    dsp::inertia<dsp::exponential_ramp> inertia_cutoff, inertia_resonance, inertia_gain;
    dsp::once_per_n timer;
    bool is_active;    
    mutable volatile int last_generation, last_calculated_generation;
    
    filter_module_with_inertia(float **ins, float **outs, float **params)
    : inertia_cutoff(dsp::exponential_ramp(128), 20)
    , inertia_resonance(dsp::exponential_ramp(128), 20)
    , inertia_gain(dsp::exponential_ramp(128), 1.0)
    , timer(128)
    , is_active(false)
    , last_generation(-1)
    , last_calculated_generation(-2)
    {}
    
    void calculate_filter()
    {
        float freq = inertia_cutoff.get_last();
        printf("freq=%g inr.cnt=%d timer.left=%d\n", freq, inertia_cutoff.count, timer.left);
        // XXXKF this is resonance of a single stage, obviously for three stages, resonant gain will be different
        float q    = inertia_resonance.get_last();
        int   mode = dsp::fastf2i_drm(*params[Metadata::par_mode]);
        //printf("freq = %f q = %f mode = %d\n", freq, q, mode);
        
        int inertia = dsp::fastf2i_drm(*params[Metadata::par_inertia]);
        if (inertia != inertia_cutoff.ramp.length()) {
            inertia_cutoff.ramp.set_length(inertia);
            inertia_resonance.ramp.set_length(inertia);
            inertia_gain.ramp.set_length(inertia);
        }
        
        FilterClass::calculate_filter(freq, q, mode, inertia_gain.get_last());
    }
    
    virtual void params_changed()
    {
        calculate_filter();
    }
    
    void on_timer()
    {
        int gen = last_generation;
        inertia_cutoff.step();
        inertia_resonance.step();
        inertia_gain.step();
        calculate_filter();
        last_calculated_generation = gen;
    }
    
    void activate()
    {
        params_changed();
        FilterClass::filter_activate();
        timer = dsp::once_per_n(FilterClass::srate / 1000);
        timer.start();
        is_active = true;
    }
    
    void set_sample_rate(uint32_t sr)
    {
        FilterClass::srate = sr;
    }

    
    void deactivate()
    {
        is_active = false;
    }

    uint32_t process(uint32_t offset, uint32_t numsamples, uint32_t inputs_mask, uint32_t outputs_mask) {
//        printf("sr=%d cutoff=%f res=%f mode=%f\n", FilterClass::srate, *params[Metadata::par_cutoff], *params[Metadata::par_resonance], *params[Metadata::par_mode]);
        uint32_t ostate = 0;
        numsamples += offset;
        while(offset < numsamples) {
            uint32_t numnow = numsamples - offset;
            // if inertia's inactive, we can calculate the whole buffer at once
            if (inertia_cutoff.active() || inertia_resonance.active() || inertia_gain.active())
                numnow = timer.get(numnow);
            
            if (outputs_mask & 1) {
                ostate |= FilterClass::process_channel(0, ins[0] + offset, outs[0] + offset, numnow, inputs_mask & 1);
            }
            if (outputs_mask & 2) {
                ostate |= FilterClass::process_channel(1, ins[1] + offset, outs[1] + offset, numnow, inputs_mask & 2);
            }
            
            if (timer.elapsed()) {
                on_timer();
            }
            offset += numnow;
        }
        return ostate;
    }
};

/// biquad filter module
class filter_audio_module: 
    public filter_module_with_inertia<dsp::biquad_filter_module, filter_metadata>, 
    public frequency_response_line_graph
{
    mutable float old_cutoff, old_resonance, old_mode;
public:    
    filter_audio_module()
    : filter_module_with_inertia<dsp::biquad_filter_module, filter_metadata>(ins, outs, params)
    {
        last_generation = 0;
        old_mode = old_resonance = old_cutoff = -1;
    }
    void params_changed()
    { 
        inertia_cutoff.set_inertia(*params[par_cutoff]);
        inertia_resonance.set_inertia(*params[par_resonance]);
        inertia_filter_module::params_changed(); 
    }
        
    bool get_graph(int index, int subindex, float *data, int points, cairo_iface *context, int *mode) const;
    int get_changed_offsets(int index, int generation, int &subindex_graph, int &subindex_dot, int &subindex_gridline) const;
};

/// Filterclavier --- MIDI controlled filter by Hans Baier
class filterclavier_audio_module: 
        public filter_module_with_inertia<dsp::biquad_filter_module, filterclavier_metadata>, 
        public frequency_response_line_graph
{        
    using audio_module<filterclavier_metadata>::ins;
    using audio_module<filterclavier_metadata>::outs;
    using audio_module<filterclavier_metadata>::params;

    const float min_gain;
    const float max_gain;
    
    int last_note;
    int last_velocity;
        
public:    
    filterclavier_audio_module();
    void params_changed();
    void activate();
    void set_sample_rate(uint32_t sr);
    void deactivate();
  
    /// MIDI control
    virtual void note_on(int channel, int note, int vel);
    virtual void note_off(int channel, int note, int vel);
    
    bool get_graph(int index, int subindex, float *data, int points, cairo_iface *context, int *mode) const;
    
private:
    void adjust_gain_according_to_filter_mode(int velocity);
};


#define MATH_E 2.718281828
class mono_audio_module:
    public audio_module<mono_metadata>
{
    typedef mono_audio_module AM;
    uint32_t srate;
    bool active;
    
    uint32_t clip_in, clip_outL, clip_outR;
    float meter_in, meter_outL, meter_outR;
    
    float * buffer;
    unsigned int pos;
    unsigned int buffer_size;
    float sign(float x) {
        if(x < 0) return -1.f;
        if(x > 0) return 1.f;
        return 0.f;
    }
    float _phase, _phase_sin_coef, _phase_cos_coef, _sc_level, _inv_atan_shape;
public:
    mono_audio_module();
    void params_changed();
    void activate();
    void set_sample_rate(uint32_t sr);
    void deactivate();
    uint32_t process(uint32_t offset, uint32_t numsamples, uint32_t inputs_mask, uint32_t outputs_mask);
};

class stereo_audio_module:
    public audio_module<stereo_metadata>
{
    typedef stereo_audio_module AM;
    float LL, LR, RL, RR;
    uint32_t srate;
    bool active;
    
    uint32_t clip_inL, clip_inR, clip_outL, clip_outR;
    float meter_inL, meter_inR, meter_outL, meter_outR, meter_phase;
    
    float * buffer;
    unsigned int pos;
    unsigned int buffer_size;
    float sign(float x) {
        if(x < 0) return -1.f;
        if(x > 0) return 1.f;
        return 0.f;
    }
    float _phase, _phase_sin_coef, _phase_cos_coef, _sc_level, _inv_atan_shape;
public:
    stereo_audio_module();
    void params_changed();
    void activate();
    void set_sample_rate(uint32_t sr);
    void deactivate();
    uint32_t process(uint32_t offset, uint32_t numsamples, uint32_t inputs_mask, uint32_t outputs_mask);
};

class analyzer_audio_module:
    public audio_module<analyzer_metadata>, public frequency_response_line_graph, public phase_graph_iface
{
    typedef analyzer_audio_module AM;
    uint32_t srate;
    bool active;
    int _accuracy;
    int _acc_old;
    int _scale_old;
    int _post_old;
    int _hold_old;
    int _smooth_old;
    uint32_t clip_L, clip_R;
    float meter_L, meter_R;
    
public:
    analyzer_audio_module();
    void params_changed();
    void activate();
    void set_sample_rate(uint32_t sr);
    void deactivate();
    uint32_t process(uint32_t offset, uint32_t numsamples, uint32_t inputs_mask, uint32_t outputs_mask);
    bool get_phase_graph(float ** _buffer, int * _length, int * _mode, bool * _use_fade, float * _fade, int * _accuracy, bool * _display) const;
    bool get_graph(int index, int subindex, float *data, int points, cairo_iface *context, int *mode) const;
    bool get_gridline(int index, int subindex, float &pos, bool &vertical, std::string &legend, cairo_iface *context) const;
    bool get_clear_all(int index) const;
    ~analyzer_audio_module();
    mutable int _mode_old;
    mutable bool _falling;
protected:
    static const int max_phase_buffer_size = 8192;
    int phase_buffer_size;
    float *phase_buffer;
    int fft_buffer_size;
    float *fft_buffer;
    int *spline_buffer;
    int plength;
    int ppos;
    int fpos;
    mutable fftwf_plan fft_plan;
    static const int max_fft_cache_size = 32768;
    static const int max_fft_buffer_size = max_fft_cache_size * 2;
    float *fft_inL, *fft_outL;
    float *fft_inR, *fft_outR;
    float *fft_smoothL, *fft_smoothR;
    float *fft_deltaL, *fft_deltaR;
    float *fft_holdL, *fft_holdR;
    float *fft_fallingL, *fft_fallingR;
    float *fft_freezeL, *fft_freezeR;
    mutable int lintrans;
    mutable int ____analyzer_phase_was_drawn_here;
    mutable int ____analyzer_sanitize;

};

};
#endif
