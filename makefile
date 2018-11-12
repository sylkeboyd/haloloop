cflags	=	-g -Wall

pcflags	=	

lflags	=	-lm -lopencv_imgproc320 -lopencv_highgui320 -lopencv_core320 -lopencv_videoio320 -lopencv_imgcodecs320

lib	=	-Lc:/opencv32/lib

inc	=	-Ic:/opencv32/install/include

oflags	=	-O2

cvfiles = 	haloloop.o

haloloop.exe	:	$(cvfiles)
	g++ $(cflags) -o haloloop.exe $(cvfiles) $(oflags) $(lib) $(lflags)

haloloop.o	:	haloloop.cpp
	g++ $(cflags) -c haloloop.cpp $(inc)

.PHONY	:	clean
clean	:
	rm $(cvfiles)
