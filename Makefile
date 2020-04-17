TARGET = fluid
OBJS = main.o display.o instance.o MACgrid.o navier.o

CXX = g++
CXXFLAGS = -std=c++11 -O2 -I/usr/include/freetype2
LDFLAGS = `pkg-config --libs --static glfw3` -lfreetype -lGLEW -lGL

.PHONY: all clean

all: $(TARGET)

clean:
	rm -f $(TARGET) $(OBJS)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<
