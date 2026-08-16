function saveProfileState(profileId) {
            const state = {};
            let container;
            switch (profileId) {
                case "simple_pvoigt": container = controls.profileSimplePVoigtContainer; break;
                case "tch_aniso": container = controls.profileTchAnisoContainer; break;
                case "split_pvoigt": container = controls.profileSplitPVoigtContainer; break;
            }
            if (container) {
                container.querySelectorAll('input[type="number"], input[type="checkbox"]').forEach(input => {
                    state[input.id] = (input.type === 'checkbox') ? input.checked : input.value;
                });
            }
            profileParamCache[profileId] = state;
        }

        function restoreProfileState(profileId) {
            const state = profileParamCache[profileId];
            if (!state || Object.keys(state).length === 0) return;
            let container;
            switch (profileId) {
                case "simple_pvoigt": container = controls.profileSimplePVoigtContainer; break;
                case "tch_aniso": container = controls.profileTchAnisoContainer; break;
                case "split_pvoigt": container = controls.profileSplitPVoigtContainer; break;
            }
            if (container) {
                container.querySelectorAll('input[type="number"], input[type="checkbox"]').forEach(input => {
                    if (state[input.id] !== undefined) {
                        input.type === 'checkbox' ? (input.checked = state[input.id]) : (input.value = state[input.id]);
                    }
                });
            }
        }        


        /*   ============================================================
         *   PROFILE PRESET LOAD / SAVE
         *   ============================================================
         *   Format: JSON. We pick JSON over XML because it parses natively
         *   in JS, is half the size, and diffs cleanly in git.
         *
         *   Preset shape:
         *     {
         *       "format":      "powder5-profile",
         *       "version":     1,
         *       "profileType": "simple_pvoigt" | "tch_aniso" | "split_pvoigt",
         *       "label":       "TCH synchrotron (ID22 default)",   // free text
         *       "parameters":  { "param-u": "0.04", "fit-u": false, ... }
         *     }
         *   The parameters block uses the DOM input IDs as keys, exactly
         *   like profileParamCache. Numeric values are stored as strings
         *   to preserve user formatting; checkboxes as booleans.
         *   Unknown keys are ignored on load (forward-compatible).
         */

        const PROFILE_PRESET_FORMAT = "powder5-profile";
        const PROFILE_PRESET_VERSION = 1;

        //   Map of "default file basename" → profile type. We probe these
        //   at startup so users can drop a Voigt_default.json (etc.) next
        //   to powder5.html and have the app pick it up automatically.
        const DEFAULT_PRESET_FILES = {
            "simple_pvoigt": "Voigt_default.json",
            "tch_aniso":     "TCH_default.json",
            "split_pvoigt":  "Split_default.json"
        };

        function getProfileContainer(profileId) {
            switch (profileId) {
                case "simple_pvoigt": return controls.profileSimplePVoigtContainer;
                case "tch_aniso":     return controls.profileTchAnisoContainer;
                case "split_pvoigt":  return controls.profileSplitPVoigtContainer;
                default: return null;
            }
        }

        function getProfileLabel(profileId) {
            return ({
                "simple_pvoigt": "Simple pVoigt",
                "tch_aniso":     "TCH (Size/Strain/Aniso)",
                "split_pvoigt":  "Split pVoigt (Asymmetric)"
            })[profileId] || profileId;
        }

        /**  Read the live DOM state for a given profile and return a
         *   preset-shaped object (suitable for JSON.stringify).
         */
        function buildProfilePreset(profileId) {
            const container = getProfileContainer(profileId);
            const parameters = {};
            if (container) {
                container.querySelectorAll('input[type="number"], input[type="checkbox"]').forEach(input => {
                    parameters[input.id] = (input.type === 'checkbox') ? input.checked : input.value;
                });
            }
            return {
                format:      PROFILE_PRESET_FORMAT,
                version:     PROFILE_PRESET_VERSION,
                profileType: profileId,
                label:       getProfileLabel(profileId),
                savedAt:     new Date().toISOString(),
                parameters
            };
        }

        /**  Validate a parsed object and apply it to the DOM. Returns
         *   the applied profile type on success, or null on failure.
         *   - Switches the profile selector if needed.
         *   - Pushes values into the inputs.
         *   - Updates profileParamCache so the cached state matches.
         *   - Re-renders the preview pattern.
         */
        function applyProfilePreset(preset, sourceLabel = "preset") {
            if (!preset || typeof preset !== 'object') {
                showToast(`Could not load ${sourceLabel}: not a JSON object.`, 'error');
                return null;
            }
            if (preset.format && preset.format !== PROFILE_PRESET_FORMAT) {
                showToast(`Could not load ${sourceLabel}: wrong format tag "${preset.format}".`, 'error');
                return null;
            }
            const profileId = preset.profileType;
            if (!["simple_pvoigt","tch_aniso","split_pvoigt"].includes(profileId)) {
                showToast(`Could not load ${sourceLabel}: unknown profileType "${profileId}".`, 'error');
                return null;
            }
            const container = getProfileContainer(profileId);
            if (!container) return null;
            const params = preset.parameters || {};

            //   Save the current profile's state before we switch away.
            saveProfileState(currentProfile);

            //   Switch the dropdown + visible container if needed.
            if (controls.profileSelect.value !== profileId) {
                controls.profileSelect.value = profileId;
                controls.profileSimplePVoigtContainer.classList.toggle('hidden', profileId !== 'simple_pvoigt');
                controls.profileTchAnisoContainer.classList.toggle('hidden',     profileId !== 'tch_aniso');
                controls.profileSplitPVoigtContainer.classList.toggle('hidden',  profileId !== 'split_pvoigt');
            }

            //   Apply preset values to the inputs of the target profile.
            let appliedCount = 0, unknownCount = 0;
            container.querySelectorAll('input[type="number"], input[type="checkbox"]').forEach(input => {
                if (Object.prototype.hasOwnProperty.call(params, input.id)) {
                    const v = params[input.id];
                    if (input.type === 'checkbox') {
                        input.checked = (v === true || v === 'true');
                    } else {
                        //   Coerce to string but keep numeric sanity.
                        const num = Number(v);
                        input.value = Number.isFinite(num) ? String(v) : input.value;
                    }
                    appliedCount++;
                }
            });
            //   Count unknown keys for diagnostics (forward-compat: ignored).
            for (const k of Object.keys(params)) {
                if (!container.querySelector(`#${CSS.escape(k)}`)) unknownCount++;
            }

            //   Sync the cache and current profile pointer.
            saveProfileState(profileId);
            currentProfile = profileId;

            //   Recompute & replot. The preset just changed every visible
            //   profile parameter, so we want the chart and any HKL-derived
            //   widths refreshed even if no experimental data is loaded yet.
            //   updatePreviewPattern() bails internally when prerequisites
            //   are missing (no data / no space group / no chart), so it's
            //   safe to call unconditionally.
            try { updateStephensAnisotropyUI(); } catch (_) {}
            try { invalidateHklCache(); } catch (_) {}     //   widths feed peak shapes
            try { updatePreviewPattern(); } catch (_) {}

            const tail = unknownCount > 0 ? ` (${unknownCount} unknown key${unknownCount===1?'':'s'} ignored)` : '';
            showToast(`Loaded ${getProfileLabel(profileId)}: ${appliedCount} parameter${appliedCount===1?'':'s'}${tail}.`, 'success');
            return profileId;
        }

        /**  Trigger a save of the current profile's parameters.
         *
         *   Uses window.showSaveFilePicker() (File System Access API) where
         *   available — that gives the user a real OS Save-As dialog, lets
         *   them pick the destination directory, and reuses it next time.
         *   Falls back to the legacy <a download> + window.prompt flow for
         *   browsers that lack the API (Firefox; iOS Safari; file:// pages
         *   in some setups).
         */
        async function saveCurrentProfilePreset() {
            const profileId = controls.profileSelect.value;
            //   Make sure the cache reflects what's on screen, then build.
            saveProfileState(profileId);
            const preset = buildProfilePreset(profileId);
            const json = JSON.stringify(preset, null, 2);

            //   Sensible default filename — matches the auto-load names so
            //   "Save → keep name" overwrites the active default.
            const stemDefault = ({
                "simple_pvoigt": "Voigt_default.json",
                "tch_aniso":     "TCH_default.json",
                "split_pvoigt":  "Split_default.json"
            })[profileId] || "profile.json";

            //   Path 1 (modern browsers): real OS Save-As dialog.
            if (typeof window.showSaveFilePicker === 'function') {
                try {
                    const handle = await window.showSaveFilePicker({
                        suggestedName: stemDefault,
                        //   Per-app picker id so the browser remembers the
                        //   last directory the user chose for presets.
                        id: 'powder5-profile-presets',
                        startIn: 'documents',
                        types: [{
                            description: 'Powder 5 profile preset',
                            accept: { 'application/json': ['.json'] }
                        }]
                    });
                    const writable = await handle.createWritable();
                    await writable.write(json);
                    await writable.close();
                    showToast(`Saved ${handle.name}`, 'success');
                    return;
                } catch (err) {
                    //   AbortError = user closed the dialog. Anything else
                    //   (permission denied, etc.) we report and fall through
                    //   to the legacy path so the save still completes.
                    if (err && err.name === 'AbortError') return;
                    console.warn('[profile-preset] showSaveFilePicker failed, falling back:', err);
                }
            }

            //   Path 2 (fallback): prompt for filename, anchor-download.
            const suggested = (window.prompt(
                "Save profile preset as (filename, .json will be appended if missing):",
                stemDefault
            ) || "").trim();
            if (!suggested) return;
            const fname = /\.json$/i.test(suggested) ? suggested : `${suggested}.json`;

            const blob = new Blob([json], { type: 'application/json' });
            const url  = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = fname;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
            showToast(`Saved ${fname}`, 'success');
        }

        /**  Open the file picker; on selection, parse and apply. */
        function loadProfilePresetFromFile() {
            const input = document.getElementById('profile-load-input');
            if (!input) return;
            input.value = '';   //  allow re-loading the same file
            input.onchange = (e) => {
                const file = e.target.files && e.target.files[0];
                if (!file) return;
                const reader = new FileReader();
                reader.onload  = (ev) => {
                    try {
                        const preset = JSON.parse(ev.target.result);
                        applyProfilePreset(preset, file.name);
                    } catch (err) {
                        showToast(`Could not parse ${file.name}: ${err.message}`, 'error');
                    }
                };
                reader.onerror = () => showToast(`Failed to read ${file.name}.`, 'error');
                reader.readAsText(file);
            };
            input.click();
        }

        /**  At startup, silently probe for Voigt_default.json / TCH_default.json
         *   / Split_default.json. If any are found, cache them in
         *   profileParamCache so that switching to that profile via the
         *   dropdown immediately restores the user's preferred defaults.
         *   The ACTIVE profile (the one selected on load) is also applied
         *   to its inputs.
         */
        async function loadDefaultProfilePresets() {
            const initialProfile = controls.profileSelect.value;
            for (const [profileId, fname] of Object.entries(DEFAULT_PRESET_FILES)) {
                try {
                    const resp = await fetch(fname, { cache: "no-store" });
                    if (!resp.ok) continue;     //  404 is fine — no default supplied
                    const text = await resp.text();
                    const preset = JSON.parse(text);
                    if (preset.format && preset.format !== PROFILE_PRESET_FORMAT) continue;
                    if (preset.profileType !== profileId) continue;

                    //   Cache it so a later switch to this profile picks it up.
                    profileParamCache[profileId] = preset.parameters || {};

                    //   If this is the currently visible profile, apply now.
                    if (profileId === initialProfile) {
                        const container = getProfileContainer(profileId);
                        if (container) {
                            const params = preset.parameters || {};
                            container.querySelectorAll('input[type="number"], input[type="checkbox"]').forEach(input => {
                                if (Object.prototype.hasOwnProperty.call(params, input.id)) {
                                    const v = params[input.id];
                                    if (input.type === 'checkbox') input.checked = (v === true || v === 'true');
                                    else                            input.value   = String(v);
                                }
                            });
                        }
                        try { updateStephensAnisotropyUI(); } catch (_) {}
                    }
                    console.log(`[profile-preset] Loaded default ${fname} for ${profileId}.`);
                } catch (err) {
                    //   fetch() rejects for file:// in some browsers — that's fine.
                    console.debug(`[profile-preset] No default ${fname}: ${err.message}`);
                }
            }
        }