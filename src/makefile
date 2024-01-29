cc = gcc
cxx = g++
wasmcc = emcc

coresrcdir = core
coreobjdir = bin
corelibs = -lm
corecflags = -Wall -ansi -pedantic

coresrc = $(wildcard $(coresrcdir)/*.c)
coreobj = $(patsubst $(coresrcdir)/%.c,$(coreobjdir)/%.o,$(coresrc))

editorinclude = -I "raylib-5.0/src" -I "imgui-1.90.1" -I "rlimgui" -I "core"
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

all: editor editor_web tests

editor: $(coreobj) $(imguiobj) $(rlimguiobj)
	$(cxx) -o editor/bin/editor $(editorinclude) $(editorflags) $(coreobj) $(imguiobj) $(rlimguiobj) editor/editor.cpp lib/libraylib_linux.a $(editorlibs)

tests: $(coreobj)
	$(cc) -o tests/trajectory/bin/tester $(editorinclude) $(coreobj) tests/trajectory/main.c -lm

# make this not dogshit later
editor_web: $(coreobj) $(imguiobj) $(rlimguiobj)
	$(wasmcc) $(rlimguiflags) -s ASYNCIFY -s USE_GLFW=3 $(editorinclude) lib/libraylib_emcc.a editor/editor.cpp $(coresrc) $(imguisrc) $(rlimguisrc) -lm -o editor/bin/wasm/arm.html

clean:
	rm -f bin/*.o || true
	rm -f editor/bin/editor || true
	rm -f editor/bin/*.o || true
	rm -f editor/bin/wasm/arm.* || true
	rm -f tests/trajectory/bin/tester || true
	rm -f tests/trajectory/csv/*.csv || true
	rm -f tests/trajectory/img/*.svg || true


$(coreobjdir)/%.o: $(coresrcdir)/%.c
	$(cc) -c $< $(corelibs) $(corecflags) -o $@

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
