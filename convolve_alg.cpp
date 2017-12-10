#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <time.h>
#include <algorithm>
#include <math.h>

using std::fstream;
using std::cout;
using std::endl;
using std::vector;
using std::sin;
using std::cos;
using std::max;

typedef struct RIFF_CHUNK {
	// RIFF chunk descriptor
	uint8_t 	RIFF[4];						// RIFF header
	uint32_t	ChunkSize;						// full chunk size
	uint8_t		WAVE[4];						// WAVE header
} riff_chk;

typedef struct FMT_CHUNK {
	// fmt subchunk
	uint8_t		fmt[4];							// FMT header
	uint32_t	Subchunk1Size;					// size of fmt chunk
	uint16_t	AudioFormat;					// audio format
	uint16_t	NumChannels;					// number of channels
	uint32_t	SamplesPerSec;					// sampling frequency (Hz)
	uint32_t	BytesPerSec;					// bytes per second
	uint16_t	blockAlign;						// 2=16bit mono, 4 = 16bit stereo
	uint16_t	bitsPerSample;					// number of bits per sample
} fmt_chk;

typedef struct DATA_CHUNK {
	uint8_t		Subchunk2ID[4];					// "data" string
	uint32_t	Subchunk2Size;					// sampled data length
} data_chk;

typedef struct WAV_HEADER {
	riff_chk RIFF_CH;
	fmt_chk FMT_CH;
	data_chk DATA_CH;
} wav_hdr;

#define INT2FLOAT 	(32768.0)
#define FLOAT2INT 	(32767.0)
#define PI         	3.141592653589793
#define TWO_PI     	(2.0 * PI)
#define SWAP(a,b)  	tempr=(a);(a)=(b);(b)=tempr

int getFileSize(FILE *input);
void print_meta(wav_hdr header);
vector<float> parseWavFile(const char wavFile[]);

bool isPow2(int x);
int nextPow2(int v);

void four1(vector<float> &data, int nn, int isign);
vector<float> complexMult(vector<float> &x, vector<float> &h);
vector<float> convolve(vector<float> &x, int N, vector<float> &h, int M, int P);

void writeLittleEndian(unsigned int word, int numBytes, FILE *wav);
void writeWav(char *filename, unsigned long numSamples, vector<float> &data, int sampleRate);

int main(int argc , char *argv[]) {

	if(argc<4) {
		printf("Missing arguments.\n");
		return 1;
	} else {
		clock_t t; // timer
		char *inputfile = argv[1];
		char *irfile = argv[2];
		char *outputfile = argv[3];
		
		printf("input file: %s\nIR file: %s\noutput file: %s\n\n",inputfile, irfile, outputfile);
		
		vector<float> x;
		vector<float> h;
		
		t = clock();
		x = parseWavFile(inputfile);
		float t0 = clock() - t;
		
		t = clock(); 
		h = parseWavFile(irfile);
		float t1 = clock() - t;
		x.shrink_to_fit();
		h.shrink_to_fit();
		
		unsigned int outputLength = x.size() + h.size() -1;

		// pad zeros
		unsigned int maxlen = outputLength;
		if(!isPow2(maxlen))
			maxlen = nextPow2(maxlen);

		while(x.size()!=maxlen)
			x.push_back(0.0);

		while(h.size()!=maxlen)
			h.push_back(0.0);

		vector<float> X;
		for(unsigned int i=0; i<x.size(); i++) {
			X.push_back(0.0);	// imaginary
		}
		for(unsigned int i=0, j=0; i<x.size(), j<x.size(); i++, j+=2) {
			X[j] = x[i];		// real
		}
		vector<float> H;
		for(unsigned int i=0; i<h.size(); i++) {
			H.push_back(0.0);	// imaginary
		}
		for(unsigned int i=0, j=0; i<h.size(), j<h.size(); i++, j+=2) {
			H[j] = h[i];		// real
		}
		X.shrink_to_fit();
		H.shrink_to_fit();
		
		t = clock();
		four1(X, X.size()/2, 1);
		float t2 = clock() - t;
		
		t= clock();
		four1(H, H.size()/2, 1);
		float t3 = clock() - t;
		
		// multiply DFTs
		vector<float> Y;
		t = clock();
		Y = complexMult(X, H);
		float t4 = clock() - t;
		
		t = clock();
		four1(Y, Y.size()/2, -1);
		float t5 = clock() - t;
		Y.shrink_to_fit();
		
		//delete X and H
		vector<float>().swap(X);
		vector<float>().swap(H);
		
		vector<float> y;
		for(unsigned int i=0; i<Y.size(); i+=2) {
			float temp = Y[i]/(float) outputLength;
			y.push_back(temp);
		}
		y.shrink_to_fit();

		// scale y[] back to ints
		for (unsigned int i=0; i<y.size(); i++) {
			y[i] = y[i] * FLOAT2INT;
		}

		// write input wav
		t = clock();
		writeWav(outputfile, y.size(), y, 44100);
		float t6 = clock() - t;
		
		float totalTime = t0 + t1 + t2 + t3 + t4 + t5 + t6;
		cout << "parseWavFile(inputfile): 			" << (float)t0/CLOCKS_PER_SEC << "s	" << t0/totalTime*100 << "%" << endl;
		cout << "parseWavFile(irfile): 				" << (float)t1/CLOCKS_PER_SEC << "s	" << t1/totalTime*100 << "%" << endl;
		cout << "four1(X, X.size()/2, 1): 			" << (float)t2/CLOCKS_PER_SEC << "s	" << t2/totalTime*100 << "%" << endl;
		cout << "four1(H, H.size()/2, 1): 			" << (float)t3/CLOCKS_PER_SEC << "s	" << t3/totalTime*100 << "%" << endl;
		cout << "complexMult(X, H): 				" << (float)t4/CLOCKS_PER_SEC << "s	" << t4/totalTime*100 << "%" << endl;
		cout << "four1(Y, Y.size()/2, -1): 			" << (float)t5/CLOCKS_PER_SEC << "s	" << t5/totalTime*100 << "%" << endl;
		cout << "writeWav(outputfile, y.size(), y, 44100): 	" << (float)t6/CLOCKS_PER_SEC << "s	" << t6/totalTime*100 << "%" << endl;
		cout << "-----" << endl;
		cout << "total time: " << (float)totalTime/CLOCKS_PER_SEC << " seconds" << endl;
		return 0;
	}
	
}

int getFileSize(FILE *input) {
	int fileSize = 0;
	fseek(input, 0, SEEK_END);
	fileSize = ftell(input);
	fseek(input, 0, SEEK_SET);
	return fileSize;
}

void print_meta(wav_hdr header) {
	// print wav file info
	cout << "RIFF header: 			\"" << header.RIFF_CH.RIFF[0] << header.RIFF_CH.RIFF[1] << header.RIFF_CH.RIFF[2] << header.RIFF_CH.RIFF[3] << "\"" << endl;
	cout << "chunk size: 			" << header.RIFF_CH.ChunkSize << endl;
	cout << "WAVE header: 			\"" << header.RIFF_CH.WAVE[0] << header.RIFF_CH.WAVE[1] << header.RIFF_CH.WAVE[2] << header.RIFF_CH.WAVE[3] << "\"" << endl;
	
	cout << "fmt: 				\"" << header.FMT_CH.fmt[0] << header.FMT_CH.fmt[1] << header.FMT_CH.fmt[2] << header.FMT_CH.fmt[3] << "\"" << endl;
	cout << "subchunk1 size: 		" << header.FMT_CH.Subchunk1Size << endl;
	// sampling rate
	cout << "sampling rate: 			" << header.FMT_CH.SamplesPerSec << endl;
	cout << "number of bits used: 		" << header.FMT_CH.bitsPerSample << endl;
	cout << "number of channels: 		" << header.FMT_CH.NumChannels << endl;
	cout << "number of bytes/second: 	" << header.FMT_CH.BytesPerSec << endl;
	cout << "audio format: 			" << header.FMT_CH.AudioFormat << endl;
	cout << "block align: 			" << header.FMT_CH.blockAlign << endl;
	
	cout << "data string: 			\"" << header.DATA_CH.Subchunk2ID[0] << header.DATA_CH.Subchunk2ID[1] << header.DATA_CH.Subchunk2ID[2] << header.DATA_CH.Subchunk2ID[3] << "\"" <<endl;
	cout << "data length: 			" << header.DATA_CH.Subchunk2Size << endl;
}

vector<float> parseWavFile(const char wavFile[]) {
	vector<float> signal;
	FILE * wav;
	wav = fopen(wavFile, "r");
	if(wav == NULL) {
		fprintf(stderr, "Error opening files: %s\n\n", wavFile);
	} else {
		int wavlength = 0;
		wavlength = getFileSize(wav);
		if(wavlength>0) {
			
			wav_hdr header;
			int riffSize = sizeof(header.RIFF_CH);
			int fmtSize = sizeof(header.FMT_CH);
			int dataSize = sizeof(header.DATA_CH);
			
			// parse wav file
			fread(&header.RIFF_CH, 1, riffSize, wav);
			fread(&header.FMT_CH, 1, fmtSize, wav);
			
			uint16_t dump[1];
			if(header.FMT_CH.Subchunk1Size>16) {
				fread(&dump, 1, header.FMT_CH.Subchunk1Size-16, wav);
			}
			fread(&header.DATA_CH, 1, dataSize, wav);
			
			//print_meta(header);
	
			int sampleSize = header.FMT_CH.bitsPerSample/8;					// bytes per sample
			int numSamples = header.DATA_CH.Subchunk2Size/sampleSize;		// how many samples in wav file
			int16_t buffer[1];

			for(int i=0; i<numSamples; i++) {
				fread(&buffer[0], sampleSize, 1, wav);
				signal.push_back(buffer[0] / INT2FLOAT);					// scale  to [-1, 1]
			}
		}
	}
	fclose(wav);
	return signal;
}

// https://stackoverflow.com/questions/108318/whats-the-simplest-way-to-test-whether-a-number-is-a-power-of-2-in-c
bool isPow2(int i) {
	if ( i <= 0 ) {
        return 0;
    }
    return !(i & (i-1));
}

// https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2
int nextPow2(int v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v++;
	return v;
}

void four1(vector<float> &data, int nn, int isign) {
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
		if (j > i) {
		    SWAP(data[j-1], data[i-1]);
		    SWAP(data[j], data[i]);
		}
		m = nn;
		while (m >= 2 && j > m) {
		    j -= m;
		    m >>= 1;
		}
		j += m;
    }
    
    mmax = 2;
    while (n > mmax) {
		istep = mmax << 1;
		theta = isign * (TWO_PI/ mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
		    for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j-1] - wi * data[j];
				tempi = wr * data[j] + wi * data[j-1];
				data[j-1] = data[i-1] - tempr;
				data[j] = data[i] - tempi;
				data[i-1] += tempr;
				data[i] += tempi;
		    }
		    wtemp = wr;
			wr += (wr*wpr) - (wi*wpi);
			wi += wi*wpr + (wtemp*wpi);
		}
		mmax = istep;
    }
}

vector<float> complexMult(vector<float> &x, vector<float> &h) {
	vector<float> y;
	float re_x, im_x, re_h, im_h;
	
	for(unsigned int i=0; i<x.size(); i+=2) {
		re_x = x[i];
		im_x = x[i+1];
		re_h = h[i];
		im_h = h[i+1];
		
		y.push_back((re_x*re_h) - (im_x*im_h));			// real
		//y.push_back((x[i]*h[i]) - (x[i+1]*h[i+1]));	// real part
		y.push_back((im_x*re_h) + (re_x*im_h));			// imaginary
		//y.push_back((x[i+1]*h[i]) + (x[i]*h[i+1]));	// imaginary part 
	}
	return y;
}

// http://www.cplusplus.com/forum/beginner/166954/
void writeLittleEndian(unsigned int word, int numBytes, FILE *wav) {
	unsigned buffer;
	while(numBytes>0) {
		buffer = word & 0xff;
		fwrite(&buffer, 1, 1, wav);
		numBytes--;
		word >>= 8;
	}
}

void writeWav(char *filename, unsigned long numSamples, vector<float> &data, int sampRate) {
	FILE *wav;
	unsigned int sampleRate;
	unsigned int numChannels;	
	unsigned int bytesPerSample;
	unsigned int byteRate;
	unsigned long i;
	
	numChannels = 1; // mono
	bytesPerSample = 2;
	
	if(sampRate<=0) sampleRate = 44100;
	else sampleRate = (unsigned int) sampRate;
	
	byteRate = sampleRate * numChannels * bytesPerSample;
	
	wav = fopen(filename, "w");
	if(wav == NULL) {
		fprintf(stderr, "Error opening files: %s\n\n", filename);
	} else {
		// begin writing
		
		// RIFF
		fwrite("RIFF", 1, 4, wav);
		writeLittleEndian(36 + bytesPerSample * numSamples * numChannels, 4, wav);	//ChunkSize
		fwrite("WAVE", 1, 4, wav);
		
		// FMT
		fwrite("fmt ", 1, 4, wav);
		writeLittleEndian(16, 4, wav);								//Subchunk1Size = 16
		writeLittleEndian(1, 2, wav);								// audio format
		writeLittleEndian(numChannels, 2, wav);						// number of channels
		writeLittleEndian(sampleRate, 4, wav); 						// sample rate
		writeLittleEndian(byteRate, 4, wav);						// byte rate
		writeLittleEndian(numChannels * bytesPerSample, 2, wav); 	// block align
		writeLittleEndian(8*bytesPerSample, 2, wav);				// bits per sample
		
		// data
		fwrite("data", 1, 4, wav);
		writeLittleEndian(bytesPerSample * numSamples * numChannels, 4, wav);	// data length
		for(i=0;i < numSamples; i++) {
			//cout << "samples left: " << numSamples-i+1 << endl;
			writeLittleEndian((unsigned int)(data[i]), bytesPerSample, wav);
		}
		fclose(wav);
	}
}
