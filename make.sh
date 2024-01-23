g++ \
	-ferror-limit=100 \
	-std=c++11 \
	-fpermissive \
	-DNO_FONT_AWESOME \
	-framework CoreVideo \
	-framework IOKit \
	-framework Cocoa \
	-framework GLUT \
	-framework OpenGL \
	-I "src/raylib-5.0/src/" \
	-I "src/imgui-1.90.1" \
	-I "src/rlimgui" \
	src/libraylib_macos.a \
	src/main.c \
	src/opw.c \
	"src/imgui-1.90.1/imgui.cpp" \
	"src/imgui-1.90.1/imgui_demo.cpp" \
	"src/imgui-1.90.1/imgui_draw.cpp" \
	"src/imgui-1.90.1/imgui_tables.cpp" \
	"src/imgui-1.90.1/imgui_widgets.cpp" \
	"src/rlimgui/rlImGui.cpp" \
	-lm -o bin/arm
