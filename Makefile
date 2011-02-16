#--------------------------------------------------------------------------
#   Copyright 2011 Taro L. Saito
# 
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#--------------------------------------------------------------------------*/
#--------------------------------------
# Silk Weaver Project
#
# Makefile
# Since: 2011/02/04
#
#--------------------------------------

INSTALL:=install
MVN:=mvn ${MVN_OPTS}
PREFIX:=${HOME}/local
INSTALL_DIR:=$(PREFIX)/genome-weaver

TARGET:=target
VERSION:=1.0-snapshot
PACKAGE:=$(TARGET)/genome-weaver-$(VERSION)

SRC:=$(shell find src)

$(TARGET)/genome-weaver-$(VERSION)-bin.tar.gz:  $(SRC)
	$(MVN) package

install: $(TARGET)/genome-weaver-$(VERSION)-bin.tar.gz
	mkdir -p "$(INSTALL_DIR)"
	mkdir -p "$(PREFIX)/bin"
	if [ -d "$(INSTALL_DIR)/current/lib" ]; then rm -rf "$(INSTALL_DIR)/current/lib"; fi
	tar xvfz $< -C "$(INSTALL_DIR)"
	ln -sfn "genome-weaver-$(VERSION)" "$(INSTALL_DIR)/current"
	ln -sf "../genome-weaver/current/bin/gw" "$(PREFIX)/bin/gw"
	ln -sf "../genome-weaver/current/bin/gw" "$(PREFIX)/bin/gw"

uninstall:
	rm -rf "$(INSTALL_DIR)/genome-weaver-$(VERSION)"

