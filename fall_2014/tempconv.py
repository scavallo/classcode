#Ben Blake
#9/2/13
#Temperature Conversions

#This program is designed to read in temperatures (in Fahrenheit) from an infile, convert them into temperatures in Celsius and Kelvin, and print them to an outfile.

infile = open("temp.txt", "r")
outfile = open("answers.txt", "w")

for line in infile:
    data = line.split()
    tempf = float(data[0])

#For every line in the infile, each set of characters is converted into a floating point number, which can then be used in mathematical calculations.
 
    tempc = (tempf - 32.0) * (0.555555)
    tempk = tempc + 273.15
    print >> outfile, '%5d' % tempf, "°F", '%5d' % tempc, "°C", '%5d' % tempk, "K"

#The temperature in Fahrenheit, Celsius, and Kelvin is printed to the outfile.  The %5d code simply writes the temperatures as integers and ensures the numbers in each column are lined up with each other.
