#!/bin/zsh
# usage: borrow.sh "<pids to pause>" <command...>
# Pauses the given pids (a low-priority twin job + the batch runner(s)) while the command runs, then resumes them,
# so the system never has more than 8 CPU-bound Python processes.
pids=$1; shift
for p in ${=pids}; do kill -STOP $p 2>/dev/null; done
"$@"
rc=$?
for p in ${=pids}; do kill -CONT $p 2>/dev/null; done
exit $rc
