.PHONY: test

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROFILE = default
PROJECT_NAME = mrbase-report-module

#################################################################################
# Rules
#################################################################################

## ==== utils ====
__utils__:

## Check health
check-health:
	utils/check-health.sh

## Format code
fmt:
	@files=$$(find . -iname "*.R" -not -path "./ldsc/*"); \
	echo $${files}; \
	for file in $${files}; do \
	  echo "    - $${file}"; \
	  Rscript -e "styler::style_file(\"$${file}\")"; \
	done
	@files=$$(find . -iname "*.Rmd" -not -path "./ldsc/*"); \
	echo $${files}; \
	for file in $${files}; do \
	  echo "    - $${file}"; \
	  Rscript -e "styler::style_file(\"$${file}\")"; \
	done

## Unit tests
test:
	@files=tests/test_*.R; \
	echo Files to be tested:; \
	for file in $${files}; do \
	  echo "    - $${file}"; \
	done; \
	for file in $${files}; do \
		echo "File: $${file}"; \
		Rscript $${file}; \
	done

#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

.PHONY: help
help:
	@echo "$$(tput bold)Params:$$(tput sgr0)"
	@echo
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}'
