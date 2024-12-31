#include <Wire.h>
#include <Adafruit_MPU6050.h>
#include <Adafruit_Sensor.h>
#include <BasicLinearAlgebra.h>

Adafruit_MPU6050 mpusensor;

//Gyro and Accelerometer Variables
const int MPU = 0x68;
float offsetAccX, offsetAccY, offsetAccZ;
float AccX, AccY, AccZ;
float rollRate; float pitchRate; float yawRate;
float offfsetRollRate = 0.0; float offfsetPitchRate = 0.0; float offfsetYawRate = 0.0;
float meaRoll, meaPitch; 
float rollAngle = 0.0; float pitchAngle = 0.0; 
float accX[6] = {0, 0, 0, 0, 0, 0};
float accY[6] = {0, 0, 0, 0, 0, 0};
float accZ[6] = {0, 0, 0, 0, 0, 0};
float gyroX[6] = {0, 0, 0, 0, 0, 0};
float gyroY[6] = {0, 0, 0, 0, 0, 0};
float gyroZ[6] = {0, 0, 0, 0, 0, 0};
float a[2] = {1.64927209, -0.70219636};
float b[3] = {0.01323107, 0.02646213, 0.01323107};
float Ts = 4; float dt = Ts/1000;
unsigned long int start = 0;
float g = 9.81;
// 2D Kalman Filter 
BLA::Matrix<5, 5>F; BLA::Matrix<5, 5>P;
BLA::Matrix<3, 3>Q; BLA::Matrix<5, 1>S;
BLA::Matrix<5, 5>I; BLA::Matrix<5, 3>K;
BLA::Matrix<3, 3>R; BLA::Matrix<3, 3>L; 
BLA::Matrix<3, 1>M; BLA::Matrix<3, 1>hx;
BLA::Matrix<3, 5>H; BLA::Matrix<5, 3>W;

float cp, sp, ct, st, tt;
float p, q, r;
float bx = 0; float by = 0; float bz = 0;


float* butterworth(float value, float arr[]) {
    arr[3] = value;
    arr[0] = a[0] * arr[1] + a[1] * arr[2] + b[0] * arr[3] + b[1] * arr[4] + b[2] * arr[5];
    arr[2] = arr[1];
    arr[1] = arr[0];
    arr[5] = arr[4];
    arr[4] = arr[3];
    return arr;
}


void gyro_signals(void){
  Wire.beginTransmission(MPU);
  Wire.write(0x1A);
  Wire.endTransmission();
  
  // Range of MPU Accelerometer
  Wire.beginTransmission(MPU);
  Wire.write(0x1C);
  Wire.write(0x10);
  Wire.endTransmission();
  
  //Getting Accelerations
  Wire.beginTransmission(MPU);
  Wire.write(0x3B);
  Wire.endTransmission();
  Wire.requestFrom(MPU, 6);
  
  int16_t AccXLSB = Wire.read()<<8|Wire.read();
  int16_t AccYLSB = Wire.read()<<8|Wire.read();
  int16_t AccZLSB = Wire.read()<<8|Wire.read();

  butterworth((float)AccXLSB/4096 - offsetAccX, accX);
  butterworth((float)AccYLSB/4096 - offsetAccY, accY);
  butterworth((float)AccZLSB/4096 + offsetAccZ, accZ);
  AccX = accX[0]; AccY = accY[0]; AccZ = accZ[0]; 

 //Range of GyroRates
  Wire.beginTransmission(MPU);
  Wire.write(0x1B);
  Wire.write(0x8);
  Wire.endTransmission();

  //Getting GyroRates
  Wire.beginTransmission(MPU);
  Wire.write(0x43);
  Wire.endTransmission();
  Wire.requestFrom(MPU, 6);
  
  int16_t GyroX = Wire.read()<<8|Wire.read();
  int16_t GyroY = Wire.read()<<8|Wire.read();
  int16_t GyroZ = Wire.read()<<8|Wire.read();  

  butterworth((float)GyroX/65.5 - offfsetRollRate, gyroX);
  butterworth((float)GyroY/65.5 - offfsetPitchRate, gyroY);
  butterworth((float)GyroZ/65.5 - offfsetYawRate, gyroZ);
  rollRate = gyroX[0]; pitchRate = gyroY[0]; yawRate = gyroZ[0]; 
  p = rollRate - bx; q = pitchRate - by; r = yawRate - bz;
  cp = cos(rollAngle/57.2957); sp = sin(rollAngle/57.2957);
  ct = cos(pitchAngle/57.2957); st = sin(pitchAngle/57.2957);
  tt = tan(pitchAngle/57.2957);
  S(0,0) += dt*(p + q*sp*tt + r*cp*tt);
  S(1,0) += dt*(q*cp - r*sp);
  rollAngle = S(0,0);
  pitchAngle = S(1,0);
  W = {dt, sp*tt*dt, cp*tt*dt, 0, cp*dt, -sp*dt, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  F = {1+(q*cp-r*sp)*tt*dt, dt*(q*sp+r*cp)/(ct*ct), -dt, -sp*tt*dt, -cp*tt*dt,
      -(q*sp+r*cp)*dt, 1, 0, -cp*dt, sp*dt,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1};
  H = {0, -ct, 0, 0, 0, ct*cp, -st*sp, 0, 0, 0, -ct*sp, -st*cp, 0, 0, 0};
  P = F*P*~F + W*Q*~W;
  L = H*P*~H + R;
  K = P*~H*Inverse(L);
  M = {AccX*g, AccY*g, AccZ*g};
  hx = {-g*st, g*ct*sp, g*ct*cp};
  S = S + K*(M - hx);
  rollAngle = S(0,0);
  pitchAngle = S(1,0);
  P = (I - K*H)*P;
}


void setup() {
  // put your setup code here, to run once:
  pinMode(13, OUTPUT);
  digitalWrite(13, LOW);
  Serial.begin(9600);
  
// Try to initialize!
  if (!mpusensor.begin()) {
    Serial.println("Failed to find MPU6050 chip");
    while (1) {
      delay(10);
    }
  }
  Serial.println("MPU6050 Found!");
  Wire.setClock(400000);
  Wire.begin();
  delay(250);
  Wire.beginTransmission(MPU);
  Wire.write(0x6B);
  Wire.write(0x00);
  Wire.endTransmission();

  digitalWrite(13, HIGH);
  Serial.println("Calibration is going on.....");
  int N = 2000; 
  float sumRollRate = 0; float sumPitchRate = 0; float sumYawRate = 0;
  float sumAccX = 0; float sumAccY = 0; float sumAccZ = 0; 
   for (int i = 1; i < N; i++){
    gyro_signals();
    sumRollRate += rollRate; sumPitchRate += pitchRate; sumYawRate += yawRate; 
    sumAccX += AccX; sumAccY += AccY; sumAccZ += AccZ; 
   }
   offfsetRollRate = sumRollRate/N; offfsetPitchRate = sumPitchRate/N; offfsetYawRate = sumYawRate/N;
   offsetAccX = sumAccX/N;
   offsetAccY = sumAccY/N;
   offsetAccZ = 1 - sumAccZ/N;
   Serial.println("Calibration done!!!");

  digitalWrite(13, LOW);
  rollAngle = 0; pitchAngle = 0;
  delay(1000);
  digitalWrite(13, HIGH);

  cp = cos(rollAngle/57.2957); sp = sin(rollAngle/57.2957);
  ct = cos(pitchAngle/57.2957); st = sin(pitchAngle/57.2957);
  tt = tan(pitchAngle/57.2957);
  p = rollRate - bx; q = pitchRate - by; r = yawRate - bz;
  F = {1+(q*cp*tt-r*sp*tt)*dt, (q*sp+r*cp)*dt/(cp*cp), -(q*sp+r*cp)*dt, 1};
  H = {0, -g*ct, g*ct*cp, -g*st*sp, -g*ct*sp, -g*st*cp};

  Q = {0.1, 0, 0, 
       0, 0.4, 0,  
       0, 0, 0.3};

  I = {1, 0, 0, 0, 0,
       0, 1, 0, 0, 0,
       0, 0, 1, 0, 0,
       0, 0, 0, 1, 0,
       0, 0, 0, 0, 1};

  R = {0.02, 0, 0,
       0, 0.02, 0,
       0, 0 , 0.02};

  P = {1, 0, 0, 0, 0,
       0, 1, 0, 0, 0,
       0, 0, 1, 0, 0,
       0, 0, 0, 1, 0,
       0, 0, 0, 0, 1};

  S = {0, 0, 0, 0, 0};
}

void loop() {
  // put your main code here, to run repeatedly:
  start = millis();
  gyro_signals();
  /*Serial.print(g*AccX);
  Serial.print("\t");
  Serial.print(g*AccY);
  Serial.print("\t");
  Serial.print(g*AccZ);
  Serial.print("\t");
  Serial.print(rollRate);
  Serial.print("\t");
  Serial.print(pitchRate);
  Serial.print("\t");
  Serial.println(yawRate);
  Serial.print("\t");*/
  Serial.print(-5);
  Serial.print("\t");
  Serial.print(rollAngle);
  Serial.print("\t");
  Serial.print(pitchAngle);
  Serial.print("\t");
  Serial.print(S(2, 0));
  Serial.print("\t");
  Serial.print(S(3, 0));
  Serial.print("\t");
  Serial.print(S(4, 0));
  Serial.print("\t");
  Serial.println(5);
  while ((millis() - start) < Ts);
}
