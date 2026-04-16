# Compilador y flags
CXX := g++
CXXFLAGS := -Wall -Wextra -O2 -std=c++17 -Iinclude
LDFLAGS := -lmatplot
AR := ar
ARFLAGS := rcs

# Directorios
SRC_DIR := src
BUILD_DIR := build
INCLUDE_DIR := include/physim/core
LIB_DIR := $(BUILD_DIR)
OBJ_DIR := $(BUILD_DIR)/obj

# Archivos fuente (asumiendo que están en src/)
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Nombre de la librería
LIB_NAME := libsimulacra.a
LIB_TARGET := $(LIB_DIR)/$(LIB_NAME)

# Archivos de cabecera
HEADERS := $(wildcard $(INCLUDE_DIR)/*.h)

# Regla principal
all: $(LIB_TARGET)

# Crear la librería estática
$(LIB_TARGET): $(OBJS) | $(LIB_DIR)
	$(AR) $(ARFLAGS) $@ $^

# Compilar objetos
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@

# Crear directorios si no existen
$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Limpiar archivos compilados
clean:
	rm -rf $(BUILD_DIR)

# Limpiar y recompilar
re: clean all

# Instalar la librería (opcional)
install: $(LIB_TARGET)
	@echo "Instalando librería..."
	# Aquí puedes agregar comandos para copiar a /usr/local/lib o similar
	# cp $(LIB_TARGET) /usr/local/lib/
	# cp -r include/physim /usr/local/include/

# Mostrar información
info:
	@echo "Fuentes: $(SRCS)"
	@echo "Objetos: $(OBJS)"
	@echo "Cabeceras: $(HEADERS)"
	@echo "Librería: $(LIB_TARGET)"

.PHONY: all clean re install info
