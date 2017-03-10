#define MaxHitsPerT0 200000
typedef uint8_t BYTE;
typedef uint16_t WORD;
typedef uint32_t DWORD;
typedef int32_t INT;

#define NMONABARS 16

// Internal ELSA data structure
typedef struct {
    double time;
    uint16_t integral[2];
    uint16_t detector_id;
} ELSA_BANK;

typedef struct {
  double tmean;
  double tleft;
  double tright;
  double qleft;
  double pleft;
  double qright;
  double pright;
  int barnum;
  bool iscomplete;
} MONABAR_HIT;

typedef struct {
  double qleft[NMONABARS];
  double qright[NMONABARS];
  double pleft[NMONABARS];
  double pright[NMONABARS];
  double tmean[NMONABARS];
  double qmean[NMONABARS];
  double tdiff[NMONABARS];
  double x[NMONABARS];
  double y[NMONABARS];
  double z[NMONABARS];
  int bar[NMONABARS];
  int mult;
  std::vector<ELSA_BANK> coincevts;
} MONA_EVENT;

// Standard UAC data structure
typedef struct {
  float     interpolation;
  uint64_t  position;
  DWORD     width;
  uint64_t  wavelet_start;
  uint64_t  wavelet_stop;
  INT       area;
  INT       peak;
  INT       detector_id;
  INT       integral[6];
}
QILS_BANK;

// Standard UAC data structure
typedef struct {
  float     interpolation;
  DWORD  position;
  DWORD     width;
  DWORD  wavelet_start;
  DWORD  wavelet_stop;
  INT       area;
  INT       peak;
  INT       detector_id;
  INT       integral[6];
}
PILS_BANK;

// DANCE-specific UAC data structure
typedef struct {
  ULong64_t position;
  DWORD    extras;
  WORD     width;
  WORD     detector_id;
  WORD     integral[2];
} CEVT_BANK;

// Anna data structure
typedef struct {
    uint64_t timestamp;  // timestamp
    uint32_t Ns;         // number of samples in waveform
    uint16_t sgate;      // short gate
    uint16_t lgate;      // long gate
    uint16_t baseline;   // baseline
    uint8_t board;       // board number
    uint8_t channel;     // channel number
} DEVT_BANK;

typedef struct{
	DWORD N;
	QILS_BANK P[MaxHitsPerT0];
        //short wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        short wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        float filtered_wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
        short something[10000]; // this is actually supported channels / supported length of PXXX bank
} test_struct_qils;

typedef struct{
	DWORD N;
	CEVT_BANK P[MaxHitsPerT0];
        short wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
        short something[10000]; // this is actually supported channels / supported length of PXXX bank
} test_struct_cevt;

typedef struct{
	DWORD N;
	PILS_BANK P[MaxHitsPerT0];
        //short wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        short wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        float filtered_wavelets[MaxHitsPerT0][1024]; // this ties to the wavelets with P
        short imported_peaks[256][16384]; // this is actually supported channels / supported length of PXXX bank
        short something[10000]; // this is actually supported channels / supported length of PXXX bank
} test_struct_pils;

#ifndef INCLUDE_TMidasBanksH
#define INCLUDE_TMidasBanksH

#include <stdint.h>

// This file defines the data structures written
// into MIDAS .mid files. They define the on-disk
// data format, they cannot be arbitrarily changed.

/// Event header

struct EventHeader_t {
  uint16_t fEventId;      ///< event id
  uint16_t fTriggerMask;  ///< event trigger mask
  uint32_t fSerialNumber; ///< event serial number
  uint32_t fTimeStamp;    ///< event timestamp in seconds
  uint32_t fDataSize;     ///< event size in bytes
};

/// Bank header

struct BankHeader_t {
  uint32_t fDataSize; 
  uint32_t fFlags; 
};

/// 16-bit data bank

struct Bank_t {
  char fName[4];      ///< bank name
  uint16_t fType;     ///< type of data (see midas.h TID_xxx)
  uint16_t fDataSize;
};

/// 32-bit data bank

struct Bank32_t {
  char fName[4];      ///< bank name
  uint32_t fType;     ///< type of data (see midas.h TID_xxx)
  uint32_t fDataSize;
};
struct Sclr_t{
  uint32_t mmtm;
  char fName[4];      ///< bank name
  uint32_t fType;     ///< type of data (see midas.h TID_xxx)
  uint32_t fDataSize;  
  uint32_t scaler[24];
};
#endif
