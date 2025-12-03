# Makefile para compilar LaTeX com temporarios no build/
# Usa XeLaTeX
# Exemplo:
# make TEX=main

# Nome do arquivo TEX passado pelo terminal
TEX ?= main

# Pasta para arquivos temporários
BUILD_DIR := build

# Cria a pasta build se não existir
$(shell mkdir -p $(BUILD_DIR))

PDF := $(TEX).pdf

all:
	@echo "Compilando $(TEX).tex..."
	-@xelatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) -aux-directory=$(BUILD_DIR) -include-directory=. $(TEX).tex
	-@if grep -q "Package biblatex" $(BUILD_DIR)/$(TEX).log 2>/dev/null; then \
		echo "Rodando biber..."; \
		biber $(BUILD_DIR)/$(TEX); \
	fi
	-@xelatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(TEX).tex
	-@xelatex -interaction=nonstopmode -output-directory=$(BUILD_DIR) $(TEX).tex
	@echo "Movendo PDF para pasta principal..."
	@cp $(BUILD_DIR)/$(PDF) ./

clean:
	rm -rf $(BUILD_DIR)/*.aux $(BUILD_DIR)/*.log $(BUILD_DIR)/*.bbl \
		$(BUILD_DIR)/*.blg $(BUILD_DIR)/*.toc $(BUILD_DIR)/*.nav \
		$(BUILD_DIR)/*.out $(BUILD_DIR)/*.snm $(BUILD_DIR)/*.synctex.gz
	rm -f $(PDF)
