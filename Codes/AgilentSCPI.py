import socket
import string
import time

from numpy import *


# Based on mclib by Thomas Schmid (http://github.com/tschmid/mclib)

class SCPI:
    PORT = 5025

    def __init__(self, host, port=PORT):
        self.host = host
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.connect((host, port))
        self.f = self.s.makefile("rb")

    def close(self):
        self.s.close()

    def sendcmd(self, cmd):
        self.s.send(cmd)

    def query(self, cmd):
        self.s.send(cmd)
        c = self.s.recv(1024)
        return c

    def reset(self):
        # reset and clear device
        self.s.send("*RST\n")
        self.s.send("*CLS\n")

    def FuncGenInit(self, setupFile, waves):

        cmd = ":MMEMory:LOAD:STATe 'INT:\\" + setupFile + "'\n"
        self.s.send(cmd)
        print(cmd)
        # self.s.send(":MMEMory:LOAD:STATe 'INT:\ELECTRO_ACOUSTIC.sta'\n")

        # output 1
        if waves[0] > 0:
            cmd = ":SOUR1:FREQ " + str(waves[0]) + "\n"
        self.s.send(cmd)
        print(cmd)
        if waves[1] > 0:
            cmd = ":SOUR1:BURSt ON\n"
            self.s.send(cmd)
            cmd = ":SOUR1:BURSt:NCYCles " + str(waves[1]) + "\n"
            self.s.send(cmd)
            print(cmd)
        else:
            cmd = ":SOUR1:BURSt OFF\n"
            self.s.send(cmd)
            print(cmd)
        self.s.send(":OUTP1 ON\n")

        # output 2
        if waves[2] > 0:
            cmd = ":SOUR2:FREQ " + str(waves[2]) + "\n"
            self.s.send(cmd)
            print(cmd)
        if waves[3] > 0:
            cmd = ":SOUR2:VOLT " + str(waves[3]) + "\n"
            self.s.send(cmd)
            print(cmd)
        phase = 90
        cmd = ":SOUR2:PHASe " + str(phase) + "\n"
        self.s.send(cmd)
        print(cmd)
        self.s.send(":OUTP2 ON\n")

    def FuncGenCh2Phase(self, phase):

        self.s.send(":OUTP2 OFF\n")
        self.s.send(":SOUR2:PHASe " + str(phase) + "\n")
        self.s.send(":OUTP2 ON\n")

    def FuncGenCh2Amplitude(self, amplitude):

        self.s.send(":OUTP2 OFF\n")
        self.s.send(":SOUR2:VOLT " + str(amplitude) + "\n")
        self.s.send(":OUTP2 ON\n")

    def OscInit(self, scales):
        # Set trigger mode.
        self.s.send(":TRIGger:MODE EDGE\n")

        # Set EDGE trigger parameters.
        #    	self.s.send(":TRIGger:EDGE:SOURCe CHANnel1\n")
        self.s.send(":TRIGger:EDGE:SOURCe EXTernal\n")
        self.s.send(":TRIGger:EDGE:LEVel 5.0\n")
        self.s.send(":TRIGger:EDGE:SLOPe POSitive\n")

        # Set time reference to left hand side of screen
        self.s.send(":TIM:REF LEFT\n")

        # Change oscilloscope settings with individual commands:
        # turn on all channels
        self.s.send(":CHANnel1:DISPlay 1\n")
        self.s.send(":CHANnel2:DISPlay 1\n")
        self.s.send(":CHANnel3:DISPlay 1\n")
        self.s.send(":CHANnel4:DISPlay 1\n")

        # Set input probe attenuation.
        self.s.send(":CHANnel1:PROBe 1\n")
        self.s.send(":CHANnel2:PROBe 1\n")
        self.s.send(":CHANnel3:PROBe 1\n")
        self.s.send(":CHANnel4:PROBe 1\n")

        # Set vertical scale and offset.
        cmd1 = ":CHANnel1:SCALe " + str(scales[0]) + "\n"
        cmd2 = ":CHANnel2:SCALe " + str(scales[1]) + "\n"
        cmd3 = ":CHANnel3:SCALe " + str(scales[2]) + "\n"
        cmd4 = ":CHANnel4:SCALe " + str(scales[3]) + "\n"
        self.s.send(cmd1)
        self.s.send(cmd2)
        self.s.send(cmd3)
        self.s.send(cmd4)
        #    	self.s.send(":CHANnel1:SCALe 5.0\n")
        #    	self.s.send(":CHANnel2:SCALe 0.002\n")
        #    	self.s.send(":CHANnel3:SCALe 0.5\n")
        #    	self.s.send(":CHANnel4:SCALe 0.001\n")

        self.s.send(":CHANnel1:OFFSet 0.0\n")
        self.s.send(":CHANnel2:OFFSet 0.0\n")
        self.s.send(":CHANnel3:OFFSet 0.0\n")
        self.s.send(":CHANnel4:OFFSet 0.0\n")

        # Set horizontal scale and offset.
        #        cmd1 = ":TIMebase:REFerence LEFT\n"
        cmd2 = ":TIMebase:SCALe " + str(scales[4]) + "\n"
        cmd3 = ":TIMebase:POSition " + str(scales[5]) + "\n"
        print(cmd1)
        print(cmd2)
        print(cmd3)
        self.s.send(cmd1)
        self.s.send(cmd2)
        self.s.send(cmd3)

        # Set the acquisition type.
        self.s.send(":ACQuire:TYPE AVERage\n")
        self.s.send(":ACQuire:COUNt 16384\n")

    #    	self.s.send(":SAVE:WAVeform:FORMat CSV\n")

    def OscRun(self):
        self.s.send(":RUN\n")

    def OscStop(self):
        self.s.send(":STOP\n")

    def OscWait(self):
        c = ''
        while c != '1':
            self.s.send("*OPC?\n")
            c = self.s.recv(1024)
            time.sleep(0.1)

    # save data to usb drive
    def OscSaveWave(self, filename):
        cmd = ":SAVE:WAVeform " + "'" + filename + "'\n"
        print(cmd)
        self.s.send(cmd)

    def MakeDir(self, SubDir):
        cmd = ":MDIRectory " + SubDir
        # cmd = ":MDIRectory "+SubDir+"'\n"
        print(cmd)
        self.s.send(cmd)

    def OscOffset3(self, scales):
        self.s.send(":CHANnel3:OFFSet 0.0\n")
        self.s.send(":CHANnel3:SCALe 0.2\n")
        self.s.send("MEASure:CLear\n")
        self.s.send("MEASure:VAVerage? DISPlay, CHANnel3\n")
        c = self.s.recv(1024)
        cmd = ":CHANnel3:OFFSet " + str(float(c)) + "\n"
        self.s.send(cmd)
        print(cmd)
        cmd = ":CHANnel3:SCALe " + str(scales[2]) + "\n"
        self.s.send(cmd)

    def OscOffset4(self, scales):
        self.s.send(":CHANnel4:OFFSet 0.0\n")
        self.s.send(":CHANnel4:SCALe 0.2\n")
        self.s.send("MEASure:CLear\n")
        self.s.send("MEASure:VAVerage? DISPlay, CHANnel4\n")
        c = self.s.recv(1024)
        cmd = ":CHANnel4:OFFSet " + str(float(c)) + "\n"
        self.s.send(cmd)
        print(cmd)
        cmd = ":CHANnel4:SCALe " + str(scales[3]) + "\n"
        self.s.send(cmd)

    def OscRead2Numpy(self, channel, npts):
        self.s.send(':WAVeform:POINts:MODE RAW\n')
        self.s.send(':WAVeform:FORMat ASCII\n')
        cmd = ":WAVeform:SOURce CHANnel" + str(channel) + "\n"
        self.s.send(cmd)

        cmd = ":WAVeform:POINts " + str(npts) + "\n"
        # print cmd
        self.s.send(cmd)

        # print "getting preamble"
        self.s.send(":WAVeform:PREamble?\n")
        try:
            c = self.s.recv(1024)
        except socket.timeout:
            print("time out getting preamble")
            p = string.split(c, ",")
            p = [y.strip('+').rstrip("\n") for y in p]
            preamble = array(map(float, p))

            # print "got preamble"

            self.s.send(':WAVeform:DATA?\n')
            skip = self.s.recv(10)
            buf = []
            self.s.settimeout(1)
        while True:
            try:
                data = self.s.recv(1024)
                buf.append(data)
            except socket.timeout:
                break

        self.s.settimeout(None)
        records = "".join(buf).split(',')
        records[len(records) - 1] = records[len(records) - 1].rstrip("\n")
        data = array(map(float, records))

        return preamble, data
