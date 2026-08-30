// status_bar.js
// ---------------------------------------------------------------------------
//  ONE PROGRESS BAR AND ONE COMMENT LINE, FOR EVERY LONG-RUNNING PROCESS.
//
//  There were four. Le Bail and Pawley drove #progress-bar in the results
//  header; charge flipping drove #cf-progress-bar plus #cf-iter-count and
//  #cf-r-factor inside its own tab; the Wyckoff search drove
//  #wyckoff-progress-bar plus #wyckoff-gen-count, #wyckoff-cc and
//  #wyckoff-status inside its own. Each was visible only while its tab was
//  open, so starting a search and switching to the plot left no indication
//  that anything was running at all.
//
//  This is a bus, not a fifth bar. The four call sites publish to it; it owns
//  the footer. The old elements are still written to by the existing code and
//  are still there -- they carry per-tab detail that belongs in the tab -- so
//  nothing had to be unpicked to add this.
//
//  Load after state.js and before refinement_controller.js.
// ---------------------------------------------------------------------------

'use strict';

/**
 * @typedef {'fit'|'cf'|'wyckoff'} StatusChannel
 */

const AppStatus = (() => {

    /** @type {Map<StatusChannel, {label:string, text:string, fraction:number, at:number}>} */
    const active = new Map();

    // Which channel the footer shows when more than one is live. In practice
    // the buttons prevent overlap, but a Wyckoff search launched from the
    // charge-flipping tab and a refinement started before it are both
    // reachable, and the footer has to pick one rather than flicker.
    const PRIORITY = ['fit', 'cf', 'wyckoff'];

    const LABELS = { fit: 'Refinement', cf: 'Charge flipping', wyckoff: 'Wyckoff search' };

    let els = null;
    let raf = 0;

    function dom() {
        if (els) return els;
        const bar = document.getElementById('app-status-bar-fill');
        if (!bar) return null;
        els = {
            root: document.getElementById('app-status-footer'),
            track: document.getElementById('app-status-bar'),
            fill: bar,
            text: document.getElementById('app-status-text')
        };
        return els;
    }

    /**
     * Coalesce redraws to one per frame.
     *
     * A charge-flipping run posts a progress message every cycle -- thousands
     * per second on the GPU path -- and each one used to write style.width
     * directly. Writing layout-affecting style from a message handler that
     * fires faster than the compositor forces synchronous style recalculation
     * and makes the UI jerkier the faster the computation goes.
     */
    function schedule() {
        if (raf) return;
        raf = requestAnimationFrame(() => { raf = 0; render(); });
    }

    function render() {
        const e = dom();
        if (!e) return;

        let chan = null;
        for (const c of PRIORITY) if (active.has(c)) { chan = c; break; }

        if (!chan) {
            e.root.classList.remove('busy');
            e.fill.style.width = '0%';
            e.text.textContent = '';
            e.text.removeAttribute('title');
            return;
        }

        const s = active.get(chan);
        e.root.classList.add('busy');

        // A notice is text only: there is no progress to show, and a full bar
        // sitting under "Stopping..." says the opposite of what is happening.
        if (s.notice) {
            e.track.classList.remove('indeterminate');
            e.fill.style.width = '0%';
            const nline = `${s.label}: ${s.text}`;
            e.text.textContent = nline;
            e.text.title = nline;
            return;
        }

        // A negative or non-finite fraction means "running, no estimate":
        // show the track striped rather than claim a percentage that would be
        // a guess. Some phases genuinely cannot report one -- building the
        // Wyckoff tables, compiling a shader -- and a bar frozen at 0% reads
        // as a hang.
        if (Number.isFinite(s.fraction) && s.fraction >= 0) {
            e.track.classList.remove('indeterminate');
            e.fill.style.width = `${Math.max(0, Math.min(100, s.fraction * 100)).toFixed(1)}%`;
        } else {
            e.track.classList.add('indeterminate');
            e.fill.style.width = '100%';
        }

        const line = s.text ? `${s.label}: ${s.text}` : `${s.label}\u2026`;
        e.text.textContent = line;
        e.text.title = line;   // the footer is narrow; the full line on hover
    }

    return {
        /**
         * Announce that a process has started.
         * @param {StatusChannel} channel
         * @param {string} [text] Initial comment.
         */
        begin(channel, text = '') {
            active.set(channel, {
                label: LABELS[channel] || channel,
                text, fraction: -1, at: Date.now()
            });
            schedule();
        },

        /**
         * Update a running process. Starts it if begin() was missed, so a
         * worker that posts progress before its start message cannot leave
         * the footer empty for the whole run.
         *
         * @param {StatusChannel} channel
         * @param {number} fraction 0..1, or a negative value for "no estimate".
         * @param {string} [text] Comment, e.g. "CC 0.9950, restart 2/3".
         */
        set(channel, fraction, text) {
            const s = active.get(channel);
            if (!s) {
                active.set(channel, {
                    label: LABELS[channel] || channel,
                    text: text || '', fraction, at: Date.now()
                });
            } else {
                s.fraction = fraction;
                if (text !== undefined) s.text = text;
            }
            schedule();
        },

        /**
         * A process has finished, been cancelled, or failed.
         * @param {StatusChannel} channel
         * @param {string} [finalText] Shown briefly before the footer clears.
         */
        end(channel, finalText) {
            if (finalText && active.has(channel)) {
                const s = active.get(channel);
                s.text = finalText;
                s.fraction = 1;
                schedule();
                setTimeout(() => { active.delete(channel); schedule(); }, 2500);
                return;
            }
            active.delete(channel);
            schedule();
        },

        /**
         * A one-shot message with no progress bar: "Stopping...", "complete --
         * 3 solutions", "the worker did not answer in time".
         *
         * These used to live in #wyckoff-status inside the Wyckoff tab. That
         * element is gone, and several of the messages had NO other channel --
         * no toast, no log line -- so without somewhere to land they would
         * simply have stopped being said.
         *
         * It does not disturb a running channel: if something is still
         * working, the notice waits rather than replacing a live progress
         * line with a stale one.
         *
         * @param {StatusChannel} channel
         * @param {string} text
         * @param {number} [ms=6000] How long to leave it up.
         */
        notice(channel, text, ms = 6000) {
            if (!text) return;
            if (active.has(channel)) { active.get(channel).text = text; schedule(); return; }
            active.set(channel, {
                label: LABELS[channel] || channel,
                text, fraction: 1, at: Date.now(), notice: true
            });
            schedule();
            setTimeout(() => {
                const s = active.get(channel);
                if (s && s.notice && s.text === text) { active.delete(channel); schedule(); }
            }, ms);
        },

        /** @returns {boolean} Whether anything is running. */
        busy() { return active.size > 0; }
    };
})();

if (typeof window !== 'undefined') window.AppStatus = AppStatus;
