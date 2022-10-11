#ifndef dsplib
#define dsplib

void BufferPushRight(int bufferSize, g_real_type *bufferArray, g_real_type newValue) {
  int outIndex;
  for (outIndex = bufferSize - 1; outIndex > 0; --outIndex) {
    bufferArray[outIndex] = bufferArray[outIndex - 1];
  }
  bufferArray[0] = newValue;
  return;
}

//фильтрация сигнала во временной области по рекурсивному соотношению
int RecursiveFilterProc(
  char reset,
  int delay,
  void*    nominatorCoefficients,
  void*    denominatorCoefficients,
  void*    inputBuffer,
  void*    outputBuffer,
  g_real_type inputValue,
  g_real_type *outputValue) {
  
  int i = 0;

  if (reset) {
    int cellIndex;
    for (cellIndex = 0; cellIndex < delay; ++cellIndex) {
	  ((g_real_type*)inputBuffer)[cellIndex] = 0;
	  ((g_real_type*)outputBuffer)[cellIndex] = 0;
    }
    return 0;
  }

  BufferPushRight(delay, (g_real_type*)inputBuffer, inputValue);
  BufferPushRight(delay, (g_real_type*)outputBuffer, 0);

  for (i = 0; i < delay; ++i) {
	((g_real_type*)outputBuffer)[0] += ((g_real_type*)inputBuffer)[i] * ((g_real_type*)nominatorCoefficients)[i];
  }
  for (i = 1; i < delay; ++i) {
	((g_real_type*)outputBuffer)[0] -= ((g_real_type*)outputBuffer)[i] * ((g_real_type*)denominatorCoefficients)[i];
  }

  *outputValue = ((g_real_type*)outputBuffer)[0];

  return 0;
}

g_real_type SosFilterStep(g_real_type* KSos, g_real_type* Buf, g_real_type X) {

  g_real_type Y;

  g_real_type v0 = Buf[0];
  g_real_type v1 = Buf[1];
  g_real_type v2 = Buf[2];

  g_real_type a0 = KSos[3];
  g_real_type a1 = KSos[4] / a0;
  g_real_type a2 = KSos[5] / a0;
  g_real_type b0 = KSos[0] / a0;
  g_real_type b1 = KSos[1] / a0;
  g_real_type b2 = KSos[2] / a0;

  v0 = X - a1 * v1 - a2 * v2;
  Y = b0 * v0 + b1 * v1 + b2 * v2;
  v2 = v1;
  v1 = v0;

  Buf[0] = v0;
  Buf[1] = v1;
  Buf[2] = v2;

  return Y;
}

g_real_type SosFilt(int SOSCount, g_real_type* Sos, g_real_type* Buf, g_real_type X) {
  g_real_type Xj;
  int j;

  for (j = 0; j < SOSCount; j++) {
	Xj = SosFilterStep(&Sos[j*6], &Buf[j*3], X);
	X = Xj;
  }

  return Xj;
}

int SOSFilterProc(
                  char reset,
                  int SOSCount,
				  void* Sos,
				  void* Buf,
				  g_real_type inputValue,
				  g_real_type *outputValue) {

  int i;

  if (reset) {
	for (i = 0; i < SOSCount; ++i) {
	  ((g_real_type*)Buf)[i*3]   = 0;
	  ((g_real_type*)Buf)[i*3+1] = 0;
	  ((g_real_type*)Buf)[i*3+2] = 0;
	}
    return 0;
  }

  *outputValue = SosFilt(SOSCount, (g_real_type*)Sos, (g_real_type*)Buf, inputValue);

  return 0;
}

#endif


