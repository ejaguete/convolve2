#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <time.h>

using std::fstream;
using std::cout;
using std::endl;
using std::vector;

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

#define INT2FLOAT (32768.0)
#define FLOAT2INT (32767.0)

int getFileSize(FILE *input);
void print_meta(wav_hdr header);
vector<float> parseWavFile(const char wavFile[]);

vector<float> convolve(vector<float> &x, int N, vector<float> &h, int M, int P);

void writeLittleEndian(unsigned int word, int numBytes, FILE *wav);
void writeWav(char *filename, unsigned long numSamples, vector<float> &data, int sampleRate);

int main(int argc , char *argv[]) {

	if(argc<4) {
		printf("Missing arguments.\n");
		return 1;
	}	else {
		char *inputfile = argv[1];
		char *irfile = argv[2];
		char *outputfile = argv[3];
		
		printf("input file: %s\nIR file: %s\noutput file: %s\n\n",inputfile, irfile, outputfile);
		
		clock_t timer;
		timer = clock();
		vector<float> x;
		vector<float> h;
		x = parseWavFile(inputfile);
		h = parseWavFile(irfile);
	
		vector<float> y;
		int P = x.size() + h.size() -1;
		y = convolve(x, x.size(), h, h.size(), P);

		// scale y[] back to ints
		for (int i=0; i<y.size(); i++) {
			y[i] = y[i] * FLOAT2INT;
		}
		
		// write input wav
		writeWav(outputfile, y.size(), y, 44100);
		timer = clock() - timer;
		printf("WAV file created in %.2f minutes\n", (float) timer/CLOCKS_PER_SEC/60);
		
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
				signal.push_back(buffer[0] / INT2FLOAT);		// scale  to [-1, 1]
			}
		}
	}
	fclose(wav);
	return signal;
}

vector<float> convolve(vector<float> &x, int N, vector<float> &h, int M, int P) {
	int n,m;
	vector<float> y;
	// clear y
	for(n=0; n<P; n++) {
		y.push_back(0.0);
	}
	
	for(n=0; n<N; n++) {
		//cout << "sample " << n+1 << " out of " << N << endl;
		for(m=0; m<M; m++) {
			y[n+m] += x[n] * h[m];
		}
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
