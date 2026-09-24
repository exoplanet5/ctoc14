#!/bin/zsh
cd /Users/mickey/solarsystem/ctoc14
while ps -eo args | grep "export_t10b.sh" | grep -v grep > /dev/null; do sleep 30; done
results/newgen/export_t10c.sh
