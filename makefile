#based on https://www.topbug.net/blog/2019/10/28/makefile-template-for-a-shared-library-in-c-with-explanations/

SHELL = /bin/sh

NAME 							:= awfmindex
MAJOR_VERSION_NUM := 0
MINOR_VERSION_NUM := 1
VERSION 					:= $(MAJOR).$(MINOR)
TARGET_LIB				:= lib$(NAME).so

BUILD_DIR 			:= build#dir to build to, deleted with clean
SRC_DIR 				:= src#source .c file dir
SO_OUTPUT_DIR		:= bin
SRCS 						:= $(wildcard $(SRC_DIR)/*.c)
# OBJECTS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRCS:.$(.c)=.o))
OBJS 						:= $(patsubst $(SRC_DIR)/%, $(BUILD_DIR)/%, $(SRCS:.c=.o))
SO_OUTPUT_FILE	:= $(SO_OUTPUT_DIR)/$(TARGET_LIB)

CC 				:= gcc
CFLAGS 		:= -std=c11 -fpic -O3 -mtune=native -mavx2 -Wall -Werror -Wextra -fopenmp -ldivsufsort64
LDFLAGS		:= -shared # linking flags


.PHONY: all
all: $(BIN_DIR) $(SO_OUTPUT_DIR) $(BUILD_DIR) $(OBJS) $(SO_OUTPUT_FILE)

#uses the object files to build the shared library into the bin directory
$(SO_OUTPUT_FILE): $(OBJS)
	$(CC) $(LDFLAGS)-o $@ $^

#create the object files from each c file in the src directory
$(BUILD_DIR)/%.o:$(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@


#create the directories if needed
$(BUILD_DIR):
	mkdir $(BUILD_DIR)

$(BIN_DIR):
	mkdir $(BIN_DIR)

$(SO_OUTPUT_DIR):
	mkdir $(SO_OUTPUT_DIR)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(BUILD_DIR) $(BIN_DIR)
	# rm $(SO_OUTPUT_FILE) $(OBJECTS)
