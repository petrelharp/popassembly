# Downloaded 10/10//2015 from https://faculty.washington.edu/browning/beagle/beagle.html#download
for x in beagle_4.1_03Oct15.pdf beagle.12Oct15.b2c.jar run.beagle.12Oct15.b2c.example release_notes
do
    wget https://faculty.washington.edu/browning/beagle/$x -O $x
done

ln -s $PWD/beagle.12Oct15.b2c.jar beagle.jar
