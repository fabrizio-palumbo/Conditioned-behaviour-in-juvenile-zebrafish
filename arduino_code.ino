bool channel1 = false;
bool channel2 = false;
bool channel3 = false;
bool channel4 = false;
bool channel5 =false ;
bool channel6 =false;
bool channelbase = false;
bool set_amplitude = false;
bool set_frequency = false;
bool set_duration = false;
int duration = 0;
int amplitude = 1;
int frequency = 1;
int deltashock= 10;//in milliseconds
int restperiod=740;// in milliseconds, delay between two shock!!!!

int pin1 = 22;
int pin2 = 26; 
int pin3 = 38;
int pin4 = 34;
int pin5 = 50;
int pin6 = 52;
int pinbase = 47; //This is the pin used for the ground! 
int clock1 = 0;
int clock2 = 0;
int clock3 = 0;
int clock4 = 0;
int clock5 = 0;
int clock6 = 0;

bool evald = false;

int bb = 0;

char  incomingByte[1024];
int serial_position = 0;
boolean stringComplete = false;

void startTimer (Tc *tc, uint32_t channel, IRQn_Type irq, uint32_t frequency) {
  pmc_set_writeprotect(false);
  pmc_enable_periph_clk((uint32_t)irq);
  TC_Configure(tc, channel, TC_CMR_WAVE | TC_CMR_WAVSEL_UP_RC | TC_CMR_TCCLKS_TIMER_CLOCK4);
  uint32_t rc = VARIANT_MCK / 128 / frequency;
  TC_SetRA(tc, channel, rc / 2);
  TC_SetRC(tc, channel, rc);
  TC_Start(tc, channel);
  tc->TC_CHANNEL[channel].TC_IER = TC_IER_CPCS;
  tc->TC_CHANNEL[channel].TC_IDR = ~TC_IER_CPCS;
  NVIC_EnableIRQ(irq);
}



void setup() {


  // put your setup code here, to run once:
  pinMode(pin1 , OUTPUT);
  pinMode(pin2, OUTPUT);
  pinMode(pin3 , OUTPUT);
  pinMode(pin4 , OUTPUT);
  pinMode(pin5 , OUTPUT);
  pinMode(pin6 , OUTPUT);
  pinMode(pinbase , OUTPUT);

  digitalWrite(pin1 , LOW);
  digitalWrite(pin2 , LOW);
  digitalWrite(pin3 , LOW);
  digitalWrite(pin4 , LOW);
  digitalWrite(pin5 , LOW);
  digitalWrite(pin6 , LOW);
    digitalWrite(pinbase ,LOW);
    channelbase=false;
  while (!Serial);

  SerialUSB.begin(57600);

  startTimer (TC1, 0, TC3_IRQn,1000);

}

void loop() {
  // put your main code here, to run repeatedly:
  while (!Serial);

  if (SerialUSB.available()) serialEvent();

  if (stringComplete) {
pinMode(13 , OUTPUT);
digitalWrite(13 , LOW);
    parseInputBuffer(incomingByte);

    memset(incomingByte, 0, sizeof(incomingByte));
    serial_position = 0;
    stringComplete = false;
  }

}

void parseInputBuffer(char* buff) {

  if (strcmp(buff, "Yo!\n") == 0) {
    SerialUSB.println("Hi!");
  }
  //else SerialUSB.println("Hi!");
  if (strcmp(buff, "1start\n") == 0) {
    channel1 = true;// true;
    channel2=true;
  }
  if (strcmp(buff, "1stop\n") == 0)  {
    channel1 = false;//false;
    channel2=false;
  }
  if (strcmp(buff, "2start\n") == 0) {
   // channel2 = true;//
  }
  if (strcmp(buff, "2stop\n") == 0) {
    //channel2 = false;//
  }
  if (strcmp(buff, "3start\n") == 0) {
    channel3 = true;//
  }
  if (strcmp(buff, "3stop\n") == 0) {
    channel3 = false;//
  }
  if (strcmp(buff, "4start\n") == 0) {
    channel4 = true;
    //channel5=true;
  }
  if (strcmp(buff, "4stop\n") == 0) {
  channel4 = false;
  //channel5=false;
  }
  if (strcmp(buff, "5start\n") == 0) {
    channel5 = true;
  }
  if (strcmp(buff, "5stop\n") == 0) {
    channel5 = false;
  }
  if (strcmp(buff, "6start\n") == 0) {
    channel6 = false;
  }
  if (strcmp(buff, "6stop\n") == 0) {

    channel6 = false;
  }
  if (strcmp(buff, "STOP\n") == 0) {
    channel1 = false;
    channel2 = false;
    channel3 = false;
    channel4 = false;
    channel5 = false;
    channel6 = false;
    if(!channelbase)
    {channelbase=true;
      }else{channelbase=false;}
      
  }



}
void TC3_Handler()
{
  TC_GetStatus(TC1, 0);
if(channelbase){
  digitalWrite(pinbase, HIGH);
  if (!channel1) {
    
    digitalWrite(pin1, HIGH);
    clock1 = 0;
  } else {
    if (clock1 < deltashock) {
      clock1++;
      digitalWrite(pin1,  LOW);
    }
    else {
      clock1++;
      digitalWrite(pin1, HIGH);

      if (clock1 >= restperiod) {
        clock1 = 0;
      }
    }

  }
  if (!channel2) {
    digitalWrite(pin2, HIGH);
    clock2 = 0;
  } else {
    if (clock2 < deltashock) {
      clock2++;

      digitalWrite(pin2, LOW);
    } else {
      clock2++;
        digitalWrite(pin2, HIGH);
  }
      if (clock2 >= restperiod) {
        evald = false;
        clock2 = 0;
      }
    }
 
 if (!channel3) {
    //REG_DACC_CHDR = 2;
    digitalWrite(pin3, HIGH);
    clock3 = 0;
  } else {
    if (clock3 < deltashock) {
      clock3++;
      digitalWrite(pin3, LOW);
    }
    else {
      clock3++;
      //REG_DACC_CHDR = 2;
      digitalWrite(pin3, HIGH);
      if (clock3 >= restperiod) {
        clock3 = 0;
      }
    }
  }
  if (!channel4) { 
    digitalWrite(pin4, HIGH);
    clock4 = 0;
  } else {
    if (clock4 < deltashock) {
      clock4++;
      digitalWrite(pin4, LOW);
    } else {
      clock4++;
      digitalWrite(pin4, HIGH);
      if (clock4 >= restperiod) {
        clock4 = 0;
      }
    }
  }
  if (!channel5) {
    digitalWrite(pin5, HIGH);
    clock5 = 0;
  } else {
    if (clock5 < deltashock) {
      clock5++;
      digitalWrite(pin5, LOW);
    } else {
      clock5++;
      digitalWrite(pin5, HIGH);
      if (clock5 >= restperiod) {
        clock5 = 0;
      }
    }
  }
  if (!channel6) { 
    digitalWrite(pin6, HIGH);
    clock6 = 0;
  } else {
    if (clock6 < deltashock) {
      clock6++;
      digitalWrite(pin6, LOW);
    } else {
      clock6++;
      digitalWrite(pin6, HIGH);
      if (clock6 >= restperiod) {
        clock6 = 0;
      }
    }
  }
}





else{digitalWrite(pinbase,LOW);  
if (!channel1) {
    
    digitalWrite(pin1, LOW);
    clock1 = 0;
  } else {
    if (clock1 < deltashock) {
      clock1++;
      digitalWrite(pin1,  HIGH);
    }
    else {
      clock1++;
      digitalWrite(pin1, LOW);

      if (clock1 >= restperiod) {
        clock1 = 0;
      }
    }

  }
  if (!channel2) {
    digitalWrite(pin2, LOW);
    clock2 = 0;
  } else {
    if (clock2 < deltashock) {
      clock2++;

      digitalWrite(pin2, HIGH);
    } else {
      clock2++;
        digitalWrite(pin2, LOW);
  }
      if (clock2 >= restperiod) {
        evald = false;
        clock2 = 0;
      }
    }
 
 if (!channel3) {
    //REG_DACC_CHDR = 2;
    digitalWrite(pin3, LOW);
    clock3 = 0;
  } else {
    if (clock3 < deltashock) {
      clock3++;
      digitalWrite(pin3, HIGH);
    }
    else {
      clock3++;
      //REG_DACC_CHDR = 2;
      digitalWrite(pin3, LOW);
      if (clock3 >= restperiod) {
        clock3 = 0;
      }
    }
  }
   if (!channel4) { 
    digitalWrite(pin4, LOW);
    clock4 = 0;
  } else {
    if (clock4 < deltashock) {
      clock4++;
      digitalWrite(pin4, HIGH);
    } else {
      clock4++;
      digitalWrite(pin4, LOW);
      if (clock4 >= restperiod) {
        clock4 = 0;
      }
    }
  }
  if (!channel5) {
    digitalWrite(pin5, LOW);
    clock5 = 0;
  } else {
    if (clock5 < deltashock) {
      clock5++;
      digitalWrite(pin5, HIGH);
    } else {
      clock5++;
      digitalWrite(pin5, LOW);
      if (clock5 >= restperiod) {
        clock5 = 0;
      }
    }
  }
   if (!channel6) { 
    digitalWrite(pin6, LOW);
    clock6 = 0;
  } else {
    if (clock6 < deltashock) {
      clock6++;
      digitalWrite(pin6, HIGH);
    } else {
      clock6++;
      digitalWrite(pin6, LOW);
      if (clock6 >= restperiod) {
        clock6 = 0;
      }
    }
  }
}
}



void serialEvent() {
  while (SerialUSB.available()) {
    char inChar = (char)SerialUSB.read();
    incomingByte[serial_position] = inChar;
    serial_position++;
    if (inChar == '\n') stringComplete = true;
  }
}

