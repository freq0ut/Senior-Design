#include <SPI.h>
#include <tsunami.h>

void setup() {
  Tsunami.begin();
  Tsunami.setOutputMode(OUTPUT_MODE_SINE);
  Tsunami.setFrequency(0, 10000.0);
}

void loop() {
  Tsunami.setAmplitude(3000.0);
  delay(1000);
  Tsunami.setAmplitude(0.0);
  delay(2000);
}
