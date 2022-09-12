#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u
set -x


read -r -d '' -a uris < <(find configs -type f -exec grep "docker://.*[\"']" {} + |
                          sed -E "s|.*(docker://.*)[\"']|\1|" |
                          sort -u && printf '\0')

echo "uris: ${uris[*]}"

for uri in "${uris[@]}"; do
    name="$(echo "$uri" | sed -E 's|docker://(.*)|\1|' | tr  -c '[:alnum:].\n' '-').img"
    singularity pull --disable-cache -F --name "containers/cache/$name" "$uri" &> /dev/null \
    && echo "Done processing $uri..." &
done

echo "Waiting for pulls to complete..."
wait
echo "DONE!"
