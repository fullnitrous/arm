cc = gcc
cxx = g++

coresrcdir = core
coreobjdir = bin
corelibs = -lm
corecflags = -Wall -ansi -pedantic

coresrc = $(wildcard $(coresrcdir)/*.c)
coreobj = $(patsubst $(coresrcdir)/%.c,$(coreobjdir)/%.o,$(coresrc))

editorinclude = -I "raylib-5.0/src" -I "imgui-1.90.1" -I "rlimgui"
editorflags = -std=c++11 -fpermissive -DNO_FONT_AWESOME
editorlibs = -lGL -lm -pthread -ldl -lrt -lX11

imguisrcdir = imgui-1.90.1
imguiobjdir = bin
imguisrc = $(wildcard $(imguisrcdir)/*.cpp)
imguiobj = $(patsubst $(imguisrcdir)/%.cpp,$(imguiobjdir)/%.o,$(imguisrc))

rlimguisrcdir = rlimgui
rlimguiobjdir = bin
rlimguiflags = -DNO_FONT_AWESOME
rlimguiinclude = -I "raylib-5.0/src" -I "imgui-1.90.1"
rlimguisrc = $(wildcard $(rlimguisrcdir)/*.cpp)
rlimguiobj = $(patsubst $(rlimguisrcdir)/%.cpp,$(rlimguiobjdir)/%.o,$(rlimguisrc))

all: $(coreobj) $(imguiobj) $(rlimguiobj)
	$(cxx) -o bin/editor $(editorinclude) $(editorflags) $(coreobj) $(imguiobj) $(rlimguiobj) main.c lib/libraylib_linux.a $(editorlibs)

clean:
	rm -f $(coreobjdir)/*.o || true
	rm -f $(coreobjdir)/editor || true

$(coreobjdir)/%.o: $(coresrcdir)/%.c
	$(cc) -c  $< $(corelibs) -o $@ $(corecflags)

$(imguiobjdir)/%.o: $(imguisrcdir)/%.cpp
	$(cxx) -c $< -o $@ 

$(rlimguiobjdir)/%.o: $(rlimguisrcdir)/%.cpp
	$(cxx) $(rlimguiinclude) $(rlimguiflags) -c $< -o $@ 





# old shit
#g++ \
#	-ferror-limit=100 \
#	-std=c++11 \
#	-fpermissive \
#	-DNO_FONT_AWESOME \
#	-framework CoreVideo \
#	-framework IOKit \
#	-framework Cocoa \
#	-framework GLUT \
#	-framework OpenGL \
#	-I "src/raylib-5.0/src/" \
#	-I "src/imgui-1.90.1" \
#	-I "src/rlimgui" \
#	src/libraylib_macos.a \
#	src/main.c \
#	src/opw.c \
#	"src/imgui-1.90.1/imgui.cpp" \
#	"src/imgui-1.90.1/imgui_demo.cpp" \
#	"src/imgui-1.90.1/imgui_draw.cpp" \
#	"src/imgui-1.90.1/imgui_tables.cpp" \
#	"src/imgui-1.90.1/imgui_widgets.cpp" \
#	"src/rlimgui/rlImGui.cpp" \
#	-lm -o bin/arm
# emcc -s ASYNCIFY -s USE_GLFW=3 -I "src/raylib-5.0/src/" src/libraylib_emcc.a src/main.c src/opw.c -lm -o bin/arm.html
